import re
import sys
import time
import warnings
import traceback
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
SCRIPT_DIR = Path(__file__).resolve().parent
SIMULATIONS_DIR = SCRIPT_DIR.parent / "simulations" / "cutoff2.5_standard_LJ"

ROOT_DIRS = [
    str(SIMULATIONS_DIR / "eps0.50_co_phase_separated"),
    str(SIMULATIONS_DIR / "eps0.25_distinct_condensates"),
    str(SIMULATIONS_DIR / "eps0.10_distinct_condensates"),
    str(SIMULATIONS_DIR / "eps0.40_touching_condensates"),
]


ROOT_LABELS = {
    ROOT_DIRS[0]: "cutoff2.5, e0.5",
    ROOT_DIRS[1]: "cutoff2.5, e0.25",
    ROOT_DIRS[2]: "cutoff2.5, e0.1",
    ROOT_DIRS[3]: "cutoff2.5, e0.4",
}

FRAME_RANGE = (200, 399)
OUTPUT_DIR = Path("cluster_comparison_output")

POLYMER_ATOM_TYPES = (1, 2)
RECEPTOR_ATOM_TYPES = (3, 4)

SPECIES_TAGS = ["A", "B", "C", "D", "receptor"]
RECEPTOR_DETAIL_SPECIES = ["C", "D"]

DAT_SUFFIX = {
    "nclusters": "nclusters.dat",
    "mean_cluster_size": "mean_cluster_size.dat",
    "largest_cluster_rg": "largest_cluster_rg.dat",
    "largest_cluster_fraction": "largest_cluster_fraction.dat",
}

SPECIES_ALT_SUFFIX = {
    "largest_cluster_rg": "rg.dat",
    "largest_cluster_fraction": "fraction.dat",
}

FRAME_COL_CANDIDATES = ["frame", "step", "timestep", "time"]
VALUE_COL_HINTS = {
    "nclusters": ["nclusters", "n_clusters", "num_clusters", "clusters"],
    "mean_cluster_size": ["mean_cluster_size", "mean_size", "avg_cluster_size", "meansize"],
    "largest_cluster_rg": ["rg", "radius_of_gyration", "largest_cluster_rg"],
    "largest_cluster_fraction": ["fraction", "frac", "largest_cluster_fraction"],
}

plt.rcParams.update({
    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "legend.fontsize": 12,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "figure.dpi": 100,
})


@dataclass
class TopologyResult:
    root_dir: str
    topology: str
    N: int
    stats: dict = field(default_factory=dict)
    polymer_fraction: tuple = None
    receptor_fraction: tuple = None


@dataclass
class RunLog:
    processed: list = field(default_factory=list)
    skipped: list = field(default_factory=list)
    missing_files: list = field(default_factory=list)


LOG = RunLog()


def discover_topology_folders(root_dir: Path):
    if not root_dir.exists():
        print(f"  [WARN] root dir does not exist: {root_dir}")
        return []
    found = []
    for p in sorted(root_dir.iterdir()):
        if not p.is_dir():
            continue
        m = re.search(r"N_?(\d+)", p.name)
        if m:
            found.append((int(m.group(1)), p))
        else:
            print(f"  [INFO] skipping non-topology folder: {p.name}")
    found.sort(key=lambda t: t[0])
    return found


def find_dat_file(topology_dir: Path, species_tag: str, quantity: str):
    expected = f"{species_tag}_{DAT_SUFFIX[quantity]}"
    candidate = topology_dir / expected
    if candidate.exists():
        return candidate

    alt_suffix = SPECIES_ALT_SUFFIX.get(quantity)
    if alt_suffix:
        alt_candidate = topology_dir / f"{species_tag}_{alt_suffix}"
        if alt_candidate.exists():
            return alt_candidate

    quantity_tokens = set(re.split(r"[_\-]", DAT_SUFFIX[quantity].replace(".dat", "")))
    species_token = species_tag.lower()
    best, best_score = None, 0
    for f in topology_dir.glob("*.dat"):
        tokens = set(re.split(r"[_\-]", f.stem.lower()))
        if species_token not in tokens:
            continue
        score = len(quantity_tokens & tokens)
        if score > best_score:
            best, best_score = f, score
    if best is not None and best_score >= 2:
        return best
    return None



def read_dat_auto(path: Path) -> pd.DataFrame:
    text = path.read_text(errors="ignore").strip().splitlines()
    if not text:
        raise ValueError("empty file")

    header_line = None
    data_lines = []
    for line in text:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            header_line = stripped.lstrip("#").strip()
            continue
        data_lines.append(stripped)

    if not data_lines:
        raise ValueError("no data rows found")

    rows = [re.split(r"[,\s]+", ln) for ln in data_lines]
    ncols = len(rows[0])
    if any(len(r) != ncols for r in rows):
        raise ValueError("inconsistent column counts")

    if header_line:
        cols = re.split(r"[,\s]+", header_line)
        if len(cols) != ncols:
            cols = [f"col{i}" for i in range(ncols)]
    else:
        if ncols == 2:
            cols = ["Frame", "Value"]
        else:
            cols = [f"col{i}" for i in range(ncols)]

    df = pd.DataFrame(rows, columns=cols)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(how="all")
    return df


def detect_column(df: pd.DataFrame, hints, fallback_index=None):
    lower_cols = {c: c.lower() for c in df.columns}
    for hint in hints:
        for c, lc in lower_cols.items():
            if hint in lc:
                return c
    if fallback_index is not None and fallback_index < len(df.columns):
        return df.columns[fallback_index]
    return None


def detect_frame_column(df: pd.DataFrame):
    col = detect_column(df, FRAME_COL_CANDIDATES, fallback_index=0)
    return col


def compute_stats(df: pd.DataFrame, value_col: str, frame_col: str,
                   frame_range):
    start, end = frame_range
    sub = df[(df[frame_col] >= start) & (df[frame_col] <= end)]
    vals = sub[value_col].to_numpy(dtype=float)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        return None
    return float(np.mean(vals)), float(np.std(vals)), int(vals.size)

def try_parse_receptor_details(topology_dir: Path, species_tag: str,
                                quantity: str, equil_frame: int):
    path = topology_dir / "receptor_receptor_details.dat"
    if not path.exists():
        return None
    try:
        df = read_dat_auto(path)
    except Exception:
        return None

    frame_col = detect_frame_column(df)
    if frame_col is None:
        return None

    hints = VALUE_COL_HINTS[quantity]
    tag_pattern = re.compile(rf"(^|[_\W]){re.escape(species_tag.lower())}([_\W]|$)")
    candidates = [c for c in df.columns
                  if tag_pattern.search(c.lower()) and any(h in c.lower() for h in hints)]
    if not candidates:
        return None
    return compute_stats(df, candidates[0], frame_col, equil_frame)


def load_quantity(topology_dir: Path, species_tag: str, quantity: str,
                   equil_frame: int):
    path = find_dat_file(topology_dir, species_tag, quantity)
    if path is None:
        if species_tag in RECEPTOR_DETAIL_SPECIES:
            fallback = try_parse_receptor_details(topology_dir, species_tag,
                                                    quantity, equil_frame)
            if fallback is not None:
                return fallback
            LOG.missing_files.append(
                f"{topology_dir} :: {species_tag}_{quantity} "
                f"(no standalone file and no matching column in "
                f"receptor_receptor_details.dat)")
        else:
            LOG.missing_files.append(f"{topology_dir} :: {species_tag}_{quantity}")
        return None
    try:
        df = read_dat_auto(path)
        frame_col = detect_frame_column(df)
        value_col = detect_column(df, VALUE_COL_HINTS[quantity], fallback_index=1)
        if frame_col is None or value_col is None or frame_col == value_col:
            raise ValueError(f"could not identify frame/value columns in {path.name}")
        return compute_stats(df, value_col, frame_col, equil_frame)
    except Exception as e:
        print(f"    [SKIP] corrupted/unreadable file {path}: {e}")
        LOG.skipped.append(str(path))
        return None


def parse_lammps_type_counts(lammps_data_path: Path):
    text = lammps_data_path.read_text(errors="ignore").splitlines()

    num_types = None
    for line in text:
        m = re.match(r"\s*(\d+)\s+atom types", line)
        if m:
            num_types = int(m.group(1))
            break
    if num_types is None:
        raise ValueError("could not find 'atom types' header in lammps.data")

    start = None
    for i, line in enumerate(text):
        if line.strip().split("#")[0].strip() == "Atoms":
            start = i + 1
            break
    if start is None:
        raise ValueError("no 'Atoms' section found in lammps.data")

    rows = []
    for line in text[start:]:
        s = line.split("#")[0].strip()
        if not s:
            if rows:
                break
            continue
        parts = s.split()
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            break

    if not rows:
        raise ValueError("Atoms section had no parsable rows")

    arr = np.array(rows)
    ncols = arr.shape[1]

    type_col, best_nunique = None, 0
    for cand in (1, 2, 3):
        if cand >= ncols:
            continue
        col = arr[:, cand]
        is_int = np.all(np.isclose(col, np.round(col)))
        in_range = col.min() >= 1 and col.max() <= num_types
        nunique = len(np.unique(col))
        if is_int and in_range and nunique > 1 and nunique > best_nunique:
            type_col, best_nunique = cand, nunique
    if type_col is None:
        raise ValueError("could not identify atom-type column in Atoms section")

    types = arr[:, type_col].astype(int)
    counts = {t: int(np.sum(types == t)) for t in range(1, num_types + 1)}
    return counts


def compute_combined_fractions(topology_dir: Path, equil_frame: int):
    lammps_data = topology_dir / "lammps.data"
    if not lammps_data.exists():
        LOG.missing_files.append(f"{topology_dir} :: lammps.data")
        return None, None

    try:
        counts = parse_lammps_type_counts(lammps_data)
    except Exception as e:
        print(f"    [SKIP] could not parse {lammps_data}: {e}")
        LOG.skipped.append(str(lammps_data))
        return None, None

    total_polymer = sum(counts.get(t, 0) for t in POLYMER_ATOM_TYPES)
    total_receptor = sum(counts.get(t, 0) for t in RECEPTOR_ATOM_TYPES)
    if total_polymer == 0 or total_receptor == 0:
        print(f"    [SKIP] zero polymer/receptor atoms in {lammps_data}")
        return None, None

    def load_frame_series(species_tag):
        path = find_dat_file(topology_dir, species_tag, "largest_cluster_fraction")
        if path is None:
            LOG.missing_files.append(f"{topology_dir} :: {species_tag}_largest_cluster_fraction")
            return None
        try:
            df = read_dat_auto(path)
            frame_col = detect_frame_column(df)
            value_col = detect_column(df, VALUE_COL_HINTS["largest_cluster_fraction"],
                                       fallback_index=1)
            return df[[frame_col, value_col]].rename(
                columns={frame_col: "Frame", value_col: "Fraction"})
        except Exception as e:
            print(f"    [SKIP] corrupted file {path}: {e}")
            LOG.skipped.append(str(path))
            return None

    a_series = load_frame_series("A")
    b_series = load_frame_series("B")
    receptor_series = load_frame_series("receptor")

    polymer_stats = None
    if a_series is not None and b_series is not None:
        merged = pd.merge(a_series, b_series, on="Frame", suffixes=("_A", "_B"))
        n_a = counts.get(POLYMER_ATOM_TYPES[0], 0)
        n_b = counts.get(POLYMER_ATOM_TYPES[1], 0)
        merged["atoms_in_cluster"] = (merged["Fraction_A"] * n_a
                                       + merged["Fraction_B"] * n_b)
        merged["polymer_fraction"] = merged["atoms_in_cluster"] / total_polymer
        polymer_stats = compute_stats(merged, "polymer_fraction", "Frame", equil_frame)
    elif a_series is not None or b_series is not None:
        only = a_series if a_series is not None else b_series
        polymer_stats = compute_stats(only, "Fraction", "Frame", equil_frame)
        print(f"    [WARN] only one polymer species fraction file found in "
              f"{topology_dir.name}; polymer fraction uses that species only")

    receptor_stats = None
    if receptor_series is not None:
        receptor_stats = compute_stats(receptor_series, "Fraction", "Frame", equil_frame)

    return polymer_stats, receptor_stats

def process_all():
    results = []
    for root_str in ROOT_DIRS:
        root_dir = Path(root_str)
        print(f"\nScanning root dir: {root_dir}")
        topo_list = discover_topology_folders(root_dir)
        if not topo_list:
            print("  [INFO] no N_* topology folders found")
            continue

        for N, topo_dir in topo_list:
            print(f"  -> {topo_dir.name} (N={N})")
            required_check = topo_dir / "lammps.data"
            if not required_check.exists():
                print(f"     [SKIP] missing lammps.data")
                LOG.skipped.append(str(topo_dir))
                continue

            res = TopologyResult(root_dir=root_str, topology=topo_dir.name, N=N)

            for quantity in ["nclusters", "mean_cluster_size",
                              "largest_cluster_rg", "largest_cluster_fraction"]:
                res.stats[quantity] = {}
                for species_tag in SPECIES_TAGS:
                    stat = load_quantity(topo_dir, species_tag, quantity, FRAME_RANGE)
                    res.stats[quantity][species_tag] = stat

            polymer_stats, receptor_stats = compute_combined_fractions(topo_dir, FRAME_RANGE)
            res.polymer_fraction = polymer_stats
            res.receptor_fraction = receptor_stats

            results.append(res)
            LOG.processed.append(str(topo_dir))

    return results


def write_fraction_csv(results, attr, filename):
    rows = []
    for r in results:
        stat = getattr(r, attr)
        if stat is None:
            continue
        mean, std, n = stat
        rows.append({
            "RootDir": r.root_dir,
            "Topology": r.topology,
            "N": r.N,
            "Mean": mean,
            "Std": std,
            "FramesUsed": n,
        })
    if not rows:
        print(f"  [WARN] no data available for {filename}, writing empty file")
        df = pd.DataFrame(columns=["RootDir", "Topology", "N", "Mean", "Std", "FramesUsed"])
    else:
        df = pd.DataFrame(rows).sort_values(["RootDir", "N"])
    out_path = OUTPUT_DIR / filename
    df.to_csv(out_path, index=False)
    print(f"  wrote {out_path}")
    return df


def write_quantity_csv(results, quantity, filename):
    rows = []
    for r in results:
        for species_tag, stat in r.stats.get(quantity, {}).items():
            if stat is None:
                continue
            mean, std, n = stat
            rows.append({
                "RootDir": r.root_dir,
                "Topology": r.topology,
                "N": r.N,
                "Species": species_tag,
                "Mean": mean,
                "Std": std,
                "FramesUsed": n,
            })
    if not rows:
        print(f"  [WARN] no data available for {filename}, writing empty file")
        df = pd.DataFrame(columns=["RootDir", "Topology", "N", "Species", "Mean", "Std", "FramesUsed"])
    else:
        df = pd.DataFrame(rows).sort_values(["RootDir", "Species", "N"])
    out_path = OUTPUT_DIR / filename
    df.to_csv(out_path, index=False)
    print(f"  wrote {out_path}")
    return df


MARKERS = ["o", "s", "^", "D", "v", "P", "X"]


def plot_fraction_figure(df, ylabel, title, out_name):
    if df.empty:
        print(f"  [WARN] no data for {out_name}, skipping plot")
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, (root, sub) in enumerate(df.groupby("RootDir")):
        sub = sub.sort_values("N")
        label = ROOT_LABELS.get(root, root)
        ax.errorbar(sub["N"], sub["Mean"], yerr=sub["Std"],
                     marker=MARKERS[i % len(MARKERS)], markersize=8,
                     capsize=4, linewidth=2, label=label)
    ax.set_xlabel("Number of receptors, N")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.4)
    ax.legend(loc="best", frameon=True)
    fig.tight_layout()
    out_path = OUTPUT_DIR / out_name
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"  wrote {out_path}")


LINESTYLES = {"A": "--", "B": ":", "receptor": "-"}


def plot_compare_quantity(df, ylabel, title, out_name):
    if df.empty:
        print(f"  [WARN] no data for {out_name}, skipping plot")
        return
    fig, ax = plt.subplots(figsize=(9, 6.5))
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, (root, sub_root) in enumerate(df.groupby("RootDir")):
        color = color_cycle[i % len(color_cycle)]
        label_base = ROOT_LABELS.get(root, root)
        for species_tag, sub in sub_root.groupby("Species"):
            sub = sub.sort_values("N")
            ax.errorbar(sub["N"], sub["Mean"], yerr=sub["Std"],
                         marker=MARKERS[i % len(MARKERS)], markersize=7,
                         capsize=3, linewidth=2, color=color,
                         linestyle=LINESTYLES.get(species_tag, "-"),
                         label=f"{label_base} ({species_tag})")
    ax.set_xlabel("Number of receptors, N")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.4)
    ax.legend(loc="best", frameon=True, fontsize=10, ncol=1)
    fig.tight_layout()
    out_path = OUTPUT_DIR / out_name
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_compare_fraction(polymer_df, receptor_df, out_name):
    fig, ax = plt.subplots(figsize=(9, 6.5))
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    roots = sorted(set(polymer_df["RootDir"]).union(receptor_df["RootDir"]))
    for i, root in enumerate(roots):
        color = color_cycle[i % len(color_cycle)]
        label_base = ROOT_LABELS.get(root, root)
        p_sub = polymer_df[polymer_df["RootDir"] == root].sort_values("N")
        r_sub = receptor_df[receptor_df["RootDir"] == root].sort_values("N")
        if not p_sub.empty:
            ax.errorbar(p_sub["N"], p_sub["Mean"], yerr=p_sub["Std"],
                         marker=MARKERS[i % len(MARKERS)], markersize=7,
                         capsize=3, linewidth=2, color=color, linestyle="--",
                         label=f"{label_base} (polymer)")
        if not r_sub.empty:
            ax.errorbar(r_sub["N"], r_sub["Mean"], yerr=r_sub["Std"],
                         marker=MARKERS[i % len(MARKERS)], markersize=7,
                         capsize=3, linewidth=2, color=color, linestyle="-",
                         label=f"{label_base} (receptor)")
    ax.set_xlabel("Number of receptors, N")
    ax.set_ylabel("Fraction in largest cluster")
    ax.set_title("Polymer vs receptor fraction in largest cluster")
    ax.grid(True, alpha=0.4)
    ax.legend(loc="best", frameon=True, fontsize=10)
    fig.tight_layout()
    out_path = OUTPUT_DIR / out_name
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_per_species_figures(quantity_df, ylabel, title_prefix, quantity_key):
    if quantity_df.empty:
        print(f"  [WARN] no data for {quantity_key} per-species plots, skipping")
        return
    for species_tag, sub in quantity_df.groupby("Species"):
        out_name = f"{species_tag}_{quantity_key}_vs_N.png"
        plot_fraction_figure(
            sub.drop(columns=["Species"]),
            ylabel,
            f"{title_prefix} — type {species_tag}",
            out_name,
        )


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    warnings.filterwarnings("ignore")

    print("=" * 70)
    print("Cluster comparison analysis (reusing existing algorithm outputs)")
    print("=" * 70)

    results = process_all()

    if not results:
        print("\nNo topologies were successfully processed. Check ROOT_DIRS "
              "and folder naming, then re-run.")
        return

    print("\nWriting CSV summaries ...")
    poly_df = write_fraction_csv(results, "polymer_fraction",
                                  "fraction_polymer_vs_N.csv")
    recep_df = write_fraction_csv(results, "receptor_fraction",
                                   "fraction_receptor_vs_N.csv")
    fraction_df = write_quantity_csv(results, "largest_cluster_fraction",
                                      "fraction_by_type_vs_N.csv")
    nclusters_df = write_quantity_csv(results, "nclusters",
                                       "nclusters_vs_N.csv")
    meansize_df = write_quantity_csv(results, "mean_cluster_size",
                                      "mean_cluster_size_vs_N.csv")
    rg_df = write_quantity_csv(results, "largest_cluster_rg",
                                "rg_vs_N.csv")

    print("\nGenerating figures ...")
    plot_fraction_figure(recep_df, "Fraction of receptors in largest cluster",
                          "Receptor fraction vs topology size",
                          "fraction_receptor_vs_N.png")
    plot_fraction_figure(poly_df, "Fraction of polymers in largest cluster",
                          "Polymer fraction vs topology size",
                          "fraction_polymer_vs_N.png")

    plot_compare_fraction(poly_df, recep_df, "COMPARE_fraction_vs_N.png")
    plot_compare_quantity(nclusters_df, "Number of clusters",
                           "Cluster count comparison",
                           "COMPARE_nclusters_vs_N.png")
    plot_compare_quantity(meansize_df, "Mean cluster size",
                           "Mean cluster size comparison",
                           "COMPARE_mean_cluster_size_vs_N.png")
    plot_compare_quantity(rg_df, "Largest cluster Rg",
                           "Radius of gyration comparison",
                           "COMPARE_rg_vs_N.png")

    plot_per_species_figures(fraction_df, "Fraction in largest cluster",
                              "Fraction in largest cluster", "fraction")
    plot_per_species_figures(nclusters_df, "Number of clusters",
                              "Number of clusters", "nclusters")
    plot_per_species_figures(meansize_df, "Mean cluster size",
                              "Mean cluster size", "mean_cluster_size")
    plot_per_species_figures(rg_df, "Largest cluster Rg",
                              "Largest cluster Rg", "rg")

    elapsed = time.time() - t0
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Processed topologies : {len(LOG.processed)}")
    print(f"Skipped items        : {len(LOG.skipped)}")
    if LOG.skipped:
        for s in LOG.skipped[:20]:
            print(f"    - {s}")
        if len(LOG.skipped) > 20:
            print(f"    ... and {len(LOG.skipped) - 20} more")
    print(f"Missing files         : {len(LOG.missing_files)}")
    if LOG.missing_files:
        for m in LOG.missing_files[:20]:
            print(f"    - {m}")
        if len(LOG.missing_files) > 20:
            print(f"    ... and {len(LOG.missing_files) - 20} more")
    print(f"Total execution time  : {elapsed:.2f} s")
    print(f"Output directory      : {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("\nFATAL ERROR:")
        traceback.print_exc()
        sys.exit(1)
