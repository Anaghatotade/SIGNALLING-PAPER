#!/usr/bin/env python3
import re
import sys
import time
import warnings
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ovito.io import import_file
from ovito.modifiers import LoadTrajectoryModifier, ClusterAnalysisModifier

SCRIPT_DIR = Path(__file__).resolve().parent
SIMULATIONS_DIR = SCRIPT_DIR.parent / "simulations"

root_dirs = [
    str(SIMULATIONS_DIR / "LJ_eps0.5_phase_separation"),
    str(SIMULATIONS_DIR / "WCA_no_phase_separation"),
]

SKIP_N_VALUES = {0, 1}

PARTICLE_TYPES = (1, 2)
PARTICLE_LABELS = {1: "Polymer", 2: "Receptor"}
START_FRAME = 200

DATA_FILENAME = "lammps.data"
TRAJ_FILENAME = "traj.lammpstrj"

OUTPUT_DIR = Path(".")
DEFAULT_CUTOFF = 2.5 

FIG_PATHS = {t: OUTPUT_DIR / f"fraction_{PARTICLE_LABELS[t].lower()}_vs_N.png" for t in PARTICLE_TYPES}
CSV_PATHS = {t: OUTPUT_DIR / f"fraction_{PARTICLE_LABELS[t].lower()}_vs_N.csv" for t in PARTICLE_TYPES}


def discover_n_folders(root_dir):
    root = Path(root_dir)
    if not root.is_dir():
        print(f"  [WARN] Root directory not found, skipping: {root_dir}")
        return []

    n_folders = []
    for entry in root.iterdir():
        if entry.is_dir():
            m = re.match(r"^N_(\d+)$", entry.name)
            if m:
                n_value = int(m.group(1))
                if n_value in SKIP_N_VALUES:
                    continue
                n_folders.append((n_value, entry))
    n_folders.sort(key=lambda t: t[0])
    return n_folders


def extract_cutoff_from_path(root_dir):
    m = re.search(r"cutoff([\d.]+)", str(root_dir))
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    warnings.warn(
        f"Could not parse an interaction cutoff from path '{root_dir}'; "
        f"using default cutoff={DEFAULT_CUTOFF}"
    )
    return DEFAULT_CUTOFF


def root_dir_label(root_dir):
    root_str = str(root_dir)
    if "LJ_eps0.5_phase_separation" in root_str:
        return "P-P = LJ, ε=0.5"
    elif "WCA_no_phase_separation" in root_str:
        return "P-P = WCA"
    parts = Path(root_dir).parts
    if len(parts) >= 2:
        return f"{parts[-2]} / {parts[-1]}"
    return str(root_dir)


def build_pipeline(data_file, traj_file):
    pipeline = import_file(str(data_file))
    traj_mod = LoadTrajectoryModifier()
    traj_mod.source.load(str(traj_file), multiple_frames=True)
    pipeline.modifiers.append(traj_mod)
    return pipeline


def verify_atom_counts(pipeline):
    n_topology = pipeline.source.data.particles.count
    data0 = pipeline.compute(0)
    n_traj = data0.particles.count
    if n_topology != n_traj:
        print(
            f"    [WARN] Atom-count mismatch: topology has {n_topology} "
            f"atoms, trajectory frame 0 has {n_traj} atoms"
        )
        return False
    return True


def largest_cluster_fractions(data, total_counts):
    cluster_ids = np.asarray(data.particles["Cluster"])
    particle_types = np.asarray(data.particles["Particle Type"])

    in_largest = cluster_ids == 1

    fractions = {}
    for t in PARTICLE_TYPES:
        type_mask = (particle_types == t)
        count_in_cluster = int(np.count_nonzero(in_largest & type_mask))
        total = total_counts[t]
        fractions[t] = count_in_cluster / total if total > 0 else np.nan

    return fractions


def process_simulation(root_dir, n_value, folder_path):
    data_file = folder_path / DATA_FILENAME
    traj_file = folder_path / TRAJ_FILENAME

    if not data_file.is_file() or not traj_file.is_file():
        print(f"  [SKIP] N={n_value}: missing files in {folder_path}")
        return "skip", None

    cutoff = extract_cutoff_from_path(root_dir)

    try:
        pipeline = build_pipeline(data_file, traj_file)
    except Exception as exc:
        print(f"  [FAIL] N={n_value}: could not open trajectory ({exc})")
        return "fail", None

    try:
        n_frames_total = pipeline.num_frames
    except Exception as exc:
        print(f"  [FAIL] N={n_value}: could not determine frame count ({exc})")
        return "fail", None

    if n_frames_total <= START_FRAME:
        print(f"  [SKIP] N={n_value}: only {n_frames_total} frames available.")
        return "skip", None

    if not verify_atom_counts(pipeline):
        print(f"  [FAIL] N={n_value}: topology/trajectory atom-count mismatch")
        return "fail", None

    try:
        data0 = pipeline.compute(0)
        types0 = np.asarray(data0.particles["Particle Type"])
    except Exception as exc:
        print(f"  [FAIL] N={n_value}: could not read particle types ({exc})")
        return "fail", None

    total_counts = {t: int(np.count_nonzero(types0 == t)) for t in PARTICLE_TYPES}

    if sum(total_counts.values()) == 0:
        print(f"  [FAIL] N={n_value}: no valid targeted particle types found")
        return "fail", None

    pipeline.modifiers.append(
        ClusterAnalysisModifier(cutoff=cutoff, sort_by_size=True)
    )

    frame_fractions = {t: [] for t in PARTICLE_TYPES}

    for frame in range(START_FRAME, n_frames_total):
        try:
            data = pipeline.compute(frame)
            fr_dict = largest_cluster_fractions(data, total_counts)
        except Exception as exc:
            print(f"    [WARN] N={n_value}, frame {frame}: analysis failed ({exc}); skipping frame")
            continue
        
        for t in PARTICLE_TYPES:
            frame_fractions[t].append(fr_dict[t])

    if len(frame_fractions[1]) == 0:
        print(f"  [FAIL] N={n_value}: no usable frames after START_FRAME")
        return "fail", None

    result = {
        "N": n_value,
        "RootDir": str(root_dir),
        "FramesUsed": int(len(frame_fractions[1])),
    }

    for t in PARTICLE_TYPES:
        arr = np.array(frame_fractions[t], dtype=float)
        label = PARTICLE_LABELS[t]
        result[f"{label}Mean"] = float(np.nanmean(arr))
        result[f"{label}Std"] = float(np.nanstd(arr))

    log_msg = f"  [OK] N={n_value}: frames={result['FramesUsed']}"
    for t in PARTICLE_TYPES:
        label = PARTICLE_LABELS[t]
        log_msg += f", {label.lower()}={result[f'{label}Mean']:.3f}±{result[f'{label}Std']:.3f}"
    print(log_msg)

    return "ok", result


def make_plot(results_by_root, mean_key, std_key, ylabel, title, out_path):
    plt.figure(figsize=(4.6, 3.5))
    colors = plt.cm.tab10.colors

    for i, (root_dir, rows) in enumerate(results_by_root.items()):
        if not rows:
            continue
        rows_sorted = sorted(rows, key=lambda r: r["N"])
        n_vals = [r["N"] for r in rows_sorted]
        means = [r[mean_key] for r in rows_sorted]
        stds = [r[std_key] for r in rows_sorted]
        color = colors[i % len(colors)]
        plt.errorbar(
            n_vals, means, yerr=stds,
            fmt="o-", capsize=3, linewidth=1.5, markersize=4,
            color=color, ecolor=color, label=root_dir_label(root_dir),
        )

    plt.xlabel("Total Receptors", fontsize=11)
    plt.ylabel(ylabel, fontsize=11)
   ## plt.title(title, fontsize=11, fontweight="bold")
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.legend(fontsize=9)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    pdf_path = out_path.with_suffix(".pdf")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Saved figure: {out_path}")
    print(f"Saved figure: {pdf_path}")


def save_csv(results_by_root, mean_key, std_key, out_path):
    rows = []
    for root_dir, entries in results_by_root.items():
        for r in sorted(entries, key=lambda x: x["N"]):
            rows.append(
                {
                    "RootDir": root_dir,
                    "N": r["N"],
                    "Mean": r[mean_key],
                    "Std": r[std_key],
                    "FramesUsed": r["FramesUsed"],
                }
            )
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(out_path, index=False)
        print(f"Saved CSV: {out_path}")


def main():
    t_start = time.time()

    n_processed = 0
    n_skipped = 0
    n_failed = 0

    results_by_root = {root: [] for root in root_dirs}

    for root_dir in root_dirs:
        print(f"\n=== Processing root directory: {root_dir} ===")
        n_folders = discover_n_folders(root_dir)
        if not n_folders:
            print("  No N_* folders found.")
            continue

        for n_value, folder_path in n_folders:
            print(f"-> N_{n_value}")
            try:
                status, result = process_simulation(root_dir, n_value, folder_path)
            except Exception:
                print(f"  [FAIL] N={n_value}: unexpected error")
                traceback.print_exc()
                n_failed += 1
                continue

            if status == "ok":
                results_by_root[root_dir].append(result)
                n_processed += 1
            elif status == "skip":
                n_skipped += 1
            else:
                n_failed += 1

    for t in PARTICLE_TYPES:
        label = PARTICLE_LABELS[t]
        make_plot(
            results_by_root, 
            mean_key=f"{label}Mean", 
            std_key=f"{label}Std",
            ylabel=rf"$\langle F_{{{label.lower()},\ cluster}} \rangle$",
            title=f"Fraction of {label} in Largest Cluster vs Total Receptors",
            out_path=FIG_PATHS[t],
        )

    for t in PARTICLE_TYPES:
        label = PARTICLE_LABELS[t]
        save_csv(results_by_root, f"{label}Mean", f"{label}Std", CSV_PATHS[t])

    runtime = time.time() - t_start
    print("\n=== Processing summary ===")
    print(f"  Simulations processed : {n_processed}")
    print(f"  Simulations skipped   : {n_skipped}")
    print(f"  Simulations failed    : {n_failed}")
    print(f"  Total runtime         : {runtime:.1f} s")


if __name__ == "__main__":
    main()
