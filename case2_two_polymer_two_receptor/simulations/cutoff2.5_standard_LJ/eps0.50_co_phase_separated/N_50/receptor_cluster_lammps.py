#!/usr/bin/env python3

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
        description="Receptor-only cluster analysis (chains 3 & 4)."
    )

    parser.add_argument("-f", "--trajectory", required=True)
    parser.add_argument("-s", "--tpr", required=True)
    parser.add_argument("-c", "--cutoff-nm", required=True, type=float)

    parser.add_argument("--atom-names", nargs="+", default=["BB", "CA"])
    parser.add_argument("--min-pairs", type=int, default=1)

    parser.add_argument("--tmin-ps", type=float, default=0.0)
    parser.add_argument("--tmax-ps", type=float, default=None)

    parser.add_argument("--out-prefix", default="receptor_")

    parser.add_argument("--backend", choices=["serial", "OpenMP", "distopia"], default="serial")

    parser.add_argument("--time-unit", choices=["ps", "ns", "us"], default="ps")
    parser.add_argument("--movie-plane", choices=["xy", "xz", "yz"], default="xy")

    parser.add_argument("--make-html-movie", action="store_true")
    parser.add_argument("--movie-out", default=None)
    parser.add_argument("--movie-stride", type=int, default=1)
    parser.add_argument("--stride", type=int, default=1)

    parser.add_argument("--verbose", action="store_true")

    return parser.parse_args()


def time_scale(unit: str) -> float:
    return {"ps": 1.0, "ns": 1e-3, "us": 1e-6}[unit]


def build_name_selection(atom_names: list[str]) -> str:
    return " or ".join(f"name {n}" for n in atom_names if n.strip())


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


def radius_of_gyration(x, m):
    m = np.asarray(m)
    if m.sum() == 0:
        return 0.0
    com = np.average(x, axis=0, weights=m)
    rg2 = np.average(np.sum((x - com) ** 2, axis=1), weights=m)
    return float(np.sqrt(rg2))


def chain_masses_and_coms(groups):
    m, c = [], []
    for g in groups:
        if hasattr(g, "masses") and g.masses is not None:
            mass = np.asarray(g.masses, float)
            m.append(mass.sum())
            c.append(g.center_of_mass())
        else:
            m.append(len(g))
            c.append(g.positions.mean(axis=0))
    return np.array(m), np.array(c)


def choose_plane(coords, plane):
    if plane == "xy":
        return coords[:, 0], coords[:, 1], "x", "y"
    if plane == "xz":
        return coords[:, 0], coords[:, 2], "x", "z"
    return coords[:, 1], coords[:, 2], "y", "z"


def plot_timeseries(x, y, xlabel, ylabel, title, out):
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def main():
    args = parse_args()

    cutoff = args.cutoff_nm * 10.0
    sel = build_name_selection(args.atom_names)

    u = mda.Universe(args.tpr, args.trajectory)

    # -----------------------------
    # LAMMPS reference (NOT used here):
    # fix wall_receptor receptor wall/lj126 zlo 32 1.0 1.0 1.1225 zhi 35 1.0 1.0 1.1225
    # -----------------------------

    chain_groups = u.atoms.split("residue")  # fallback grouping assumption

    # FIXED: receptor chains 3 and 4 (0-based indexing)
    receptor_ids = [2, 3]
    receptor_groups = [chain_groups[i] for i in receptor_ids]

    selected = [g.select_atoms(sel) for g in receptor_groups]
    n_chains = len(selected)

    cutoff = args.cutoff_nm * 10.0
    tscale = time_scale(args.time_unit)

    times, nclusters, rg, frac, mean_size = [], [], [], [], []

    for ts in tqdm(u.trajectory[::args.stride], disable=not args.verbose):

        if ts.time < args.tmin_ps:
            continue
        if args.tmax_ps and ts.time > args.tmax_ps:
            break

        adj = build_adjacency(selected, ts.dimensions, cutoff, args.min_pairs, args.backend)
        comps = connected_components(adj)

        sizes = np.array([len(c) for c in comps])
        if len(sizes) == 0:
            continue

        largest = comps[np.argmax(sizes)]

        masses, coms = chain_masses_and_coms(selected)
        rg_val = radius_of_gyration(coms[np.array(largest)], masses[np.array(largest)])

        times.append(ts.time)
        nclusters.append(len(comps))
        rg.append(rg_val / 10.0)
        frac.append(100 * len(largest) / n_chains)
        mean_size.append(sizes.mean())

    times = np.array(times)

    prefix = args.out_prefix

    np.savetxt(prefix + "nclusters.dat", np.c_[times, nclusters])
    np.savetxt(prefix + "rg.dat", np.c_[times, rg])
    np.savetxt(prefix + "fraction.dat", np.c_[times, frac])
    np.savetxt(prefix + "mean_size.dat", np.c_[times, mean_size])

    plot_timeseries(times, nclusters, "time", "clusters", "Clusters", prefix + "nclusters.png")
    plot_timeseries(times, rg, "time", "Rg", "Rg", prefix + "rg.png")
    plot_timeseries(times, frac, "time", "% in largest cluster", "Fraction", prefix + "fraction.png")
    plot_timeseries(times, mean_size, "time", "size", "Mean size", prefix + "mean_size.png")


if __name__ == "__main__":
    main()
