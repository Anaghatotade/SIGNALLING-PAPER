from __future__ import annotations

import argparse
from collections import Counter, deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import plotly.colors as pc
import plotly.graph_objects as go
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Distance-based cluster analysis for multi-chain receptor from LAMMPS trajectory."
    )

    parser.add_argument("--lammps-data", default="lammps.data")
    parser.add_argument("--trajectory", default="traj.lammpstrj")

    parser.add_argument("-c", "--cutoff", required=True, type=float)

    parser.add_argument(
        "--atom-types",
        nargs="+",
        type=int,
        default=[3, 4],
        help="Receptor/Ligand atom types (default: 3 4)"
    )

    parser.add_argument("--min-pairs", type=int, default=1)

    parser.add_argument(
        "--chain-by",
        choices=["auto", "molnum", "resid", "resid-block", "atom-block"],
        default="auto",
    )

    parser.add_argument("--atoms-per-chain", type=int, default=None)
    parser.add_argument("--residues-per-chain", type=int, default=None)

    parser.add_argument("--stride", type=int, default=1)

    parser.add_argument(
        "--frame-mode",
        choices=["time", "index"],
        default="time"
    )

    parser.add_argument("--fmin", type=float, default=1e7)
    parser.add_argument("--fmax", type=float, default=2e7)

    parser.add_argument("--out-prefix", default="D_")

    parser.add_argument("--backend", choices=["serial", "OpenMP", "distopia"], default="serial")

    parser.add_argument("--make-html-movie", action="store_true")
    parser.add_argument("--movie-stride", type=int, default=1)
    parser.add_argument("--movie-plane", choices=["xy", "xz", "yz"], default="xy")
    parser.add_argument("--movie-out", default=None)

    parser.add_argument("--verbose", action="store_true")

    return parser.parse_args()


def build_type_selection(atom_types):
    return " or ".join(f"type {t}" for t in atom_types)


def get_chain_groups(u, chain_by="auto", residues_per_chain=None, atoms_per_chain=None):
    atoms = u.atoms

    if chain_by == "auto":
        if hasattr(atoms, "molnums"):
            chain_by = "molnum"
        elif hasattr(atoms, "resids"):
            chain_by = "resid"
        else:
            raise RuntimeError("Use explicit chain-by atom-block")

    if chain_by == "molnum":
        molnums = np.asarray(atoms.molnums)
        return [atoms[molnums == m] for m in np.unique(molnums)]

    if chain_by == "resid":
        resids = np.asarray(atoms.resids)
        return [atoms[resids == r] for r in np.unique(resids)]

    if chain_by == "resid-block":
        residues = u.residues
        n = len(residues)
        out = []
        for i in range(0, n, residues_per_chain):
            out.append(residues[i:i + residues_per_chain].atoms)
        return out

    if chain_by == "atom-block":
        n = len(atoms)
        out = []
        for i in range(0, n, atoms_per_chain):
            out.append(atoms[i:i + atoms_per_chain])
        return out

    raise RuntimeError("Invalid chain-by")


def apply_z_slab(chain_groups, zlo=32.0, zhi=35.0):
    filtered = []
    for ag in chain_groups:
        pos = ag.positions
        mask = (pos[:, 2] >= zlo) & (pos[:, 2] <= zhi)
        sel = ag[mask]
        if len(sel) > 0:
            filtered.append(sel)
    return filtered


def connected_components(adj):
    n = len(adj)
    seen = np.zeros(n, dtype=bool)
    comps = []

    for i in range(n):
        if seen[i]:
            continue
        q = deque([i])
        seen[i] = True
        comp = []
        while q:
            u = q.popleft()
            comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    q.append(v)
        comps.append(comp)

    return comps


def build_adjacency(groups, box, cutoff, min_pairs, backend):
    n = len(groups)
    adj = [[] for _ in range(n)]

    for i in range(n - 1):
        pi = groups[i].positions
        if len(pi) == 0:
            continue

        for j in range(i + 1, n):
            pj = groups[j].positions
            if len(pj) == 0:
                continue

            d = distance_array(pi, pj, box=box, backend=backend)
            if np.count_nonzero(d <= cutoff) >= min_pairs:
                adj[i].append(j)
                adj[j].append(i)

    return adj


def chain_masses_and_coms(groups):
    masses = np.zeros(len(groups))
    coms = np.zeros((len(groups), 3))

    for i, ag in enumerate(groups):
        try:
            m = np.asarray(ag.masses)
            if m.sum() > 0:
                masses[i] = m.sum()
                coms[i] = ag.center_of_mass()
                continue
        except:
            pass

        masses[i] = len(ag)
        coms[i] = ag.positions.mean(axis=0)

    return masses, coms


def minimum_image(r1, r2, box):
    dr = r2 - r1
    L = np.asarray(box[:3], float)
    dr -= L * np.round(dr / L)
    return dr


def radius_of_gyration(points, weights):
    wsum = weights.sum()
    if wsum == 0:
        return 0.0
    c = (points * weights[:, None]).sum(axis=0) / wsum
    return np.sqrt(((weights[:, None] * (points - c) ** 2).sum()) / wsum)


def unwrap_cluster(coms, cluster, box, adj):
    cluster = set(cluster)
    root = list(cluster)[0]

    placed = {root: coms[root].copy()}
    q = deque([root])

    while q:
        i = q.popleft()
        for j in adj[i]:
            if j in cluster and j not in placed:
                placed[j] = placed[i] + minimum_image(coms[i], coms[j], box)
                q.append(j)

    return np.array([placed[i] for i in cluster])


def main():
    args = parse_args()

    u = mda.Universe(args.lammps_data, args.trajectory, format="LAMMPSDUMP")

    selection = build_type_selection(args.atom_types)

    chain_groups = get_chain_groups(
        u,
        chain_by=args.chain_by,
        residues_per_chain=args.residues_per_chain,
        atoms_per_chain=args.atoms_per_chain
    )

    selected = []
    for c in chain_groups:
        s = c.select_atoms(selection)
        if len(s) > 0:
            selected.append(s)

    
    n_chains = len(selected)
    if n_chains == 0:
        raise RuntimeError("No chains after filtering")

    times = []
    nclusters_series = []
    rg_series = []
    frac_series = []

    frame_id = 0

    traj = u.trajectory[::args.stride]
    iterator = tqdm(traj, disable=not args.verbose)

    for ts in iterator:
        frame_id += 1
        t = float(ts.time)

        if args.frame_mode == "time":
            if t < args.fmin:
                continue
            if args.fmax is not None and t > args.fmax:
                break
        else:
            if frame_id < int(args.fmin):
                continue
            if frame_id > int(args.fmax):
                break

        adj = build_adjacency(selected, ts.dimensions, args.cutoff, args.min_pairs, args.backend)
        comps = connected_components(adj)

        sizes = np.array([len(c) for c in comps])

        masses, coms = chain_masses_and_coms(selected)

        largest = comps[int(np.argmax(sizes))]
        frac = 100 * len(largest) / n_chains

        coords = unwrap_cluster(coms, largest, ts.dimensions, adj)
        rg = radius_of_gyration(coords, masses[np.array(largest)])

        times.append(t)
        nclusters_series.append(len(comps))
        rg_series.append(rg)
        frac_series.append(frac)

        if args.verbose:
            iterator.set_postfix(clusters=len(comps), rg=rg, frac=frac)

    np.savetxt(args.out_prefix + "nclusters.dat", np.c_[times, nclusters_series])
    np.savetxt(args.out_prefix + "rg.dat", np.c_[times, rg_series])
    np.savetxt(args.out_prefix + "fraction.dat", np.c_[times, frac_series])

    plt.plot(times, nclusters_series)
    plt.savefig(args.out_prefix + "nclusters.png", dpi=300)
    plt.close()

    plt.plot(times, rg_series)
    plt.savefig(args.out_prefix + "rg.png", dpi=300)
    plt.close()

    plt.plot(times, frac_series)
    plt.savefig(args.out_prefix + "fraction.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
