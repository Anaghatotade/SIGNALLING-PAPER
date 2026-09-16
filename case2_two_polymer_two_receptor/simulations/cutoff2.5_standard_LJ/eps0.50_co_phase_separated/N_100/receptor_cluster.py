#!/usr/bin/env python3
"""
Distance-based cluster analysis for LAMMPS multi-chain polymer simulations.
Slab geometry with two receptors (atom types 3 & 4) near the top wall.

System layout (from diagnostic):
  type 1 : 15000 atoms  ─┐ polymer chains (100 chains × 300 beads)
  type 2 : 15000 atoms  ─┘
  type 3 :   100 atoms     receptor A  (chain index 3 by convention)
  type 4 :   100 atoms     receptor B  (chain index 4 by convention)
  Total  : 30200 atoms

Wall fix (mirrored as constants):
  fix wall_lig receptor wall/lj126 zlo 32 1.0 1.0 1.1225
                                   zhi 35 1.0 1.0 1.1225

Chain-grouping strategy
  --chain-by molnum      : uses mol column in dump (preferred if available)
  --chain-by atomtype    : groups every unique atom-type as one "chain"
                           → types 1&2 each become one super-group
                           (only useful for whole-type analysis)
  --chain-by type-molecule: NEW default — splits type-1+2 atoms into
                           equal-sized polymer chains AND keeps type-3/4
                           as individual receptor chains.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import plotly.colors as pc
import plotly.graph_objects as go
from tqdm import tqdm


# ── LAMMPS wall / receptor constants ─────────────────────────────────────────
WALL_ZLO     = 32.0    # Å  (zlo in fix wall_lig)
WALL_ZHI     = 35.0    # Å  (zhi in fix wall_lig)
WALL_EPSILON = 1.0
WALL_SIGMA   = 1.0
WALL_CUTOFF  = 1.1225  # Å  = 2^(1/6)*sigma, purely repulsive WCA

# Atom types that are receptors (LAMMPS type strings)
RECEPTOR_TYPES    = {"3", "4"}
# Atom types that make up polymer chains
POLYMER_TYPES     = {"1", "2"}
# How many atoms per polymer chain (15000 type-1 + 15000 type-2) / 100 chains
POLYMER_CHAIN_LEN = 300   # override with --polymer-chain-len


# ═════════════════════════════════════════════════════════════════════════════
#  CLI
# ═════════════════════════════════════════════════════════════════════════════

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Cluster analysis for LAMMPS slab simulations with receptors "
            "(types 3 & 4) near the top wall (zlo=32 Å, zhi=35 Å)."
        )
    )
    p.add_argument("-f", "--trajectory", required=True,
                   help="LAMMPS dump trajectory (.lammpstrj / .dump)")
    p.add_argument("-s", "--topology", required=True,
                   help="LAMMPS data file (.data)")
    p.add_argument("-c", "--cutoff-nm", required=True, type=float,
                   help="Distance cutoff in nm to connect two chains")

    p.add_argument("--atom-names", nargs="+", default=None,
                   help="Atom names for clustering. Omit to use all atoms in each chain.")
    p.add_argument("--min-pairs", type=int, default=1,
                   help="Min interchain atom pairs within cutoff to connect chains (default: 1)")

    p.add_argument("--tmin-ps", type=float, default=0.0,
                   help="Start time in ps (default: 0)")
    p.add_argument("--tmax-ps", type=float, default=None,
                   help="End time in ps (default: end)")

    p.add_argument("--out-prefix", default="receptor_",
                   help="Output file prefix (default: receptor_)")

    p.add_argument("--backend",
                   choices=["serial", "OpenMP", "distopia"], default="serial")

    p.add_argument("--chain-by",
                   choices=["molnum", "atomtype", "type-molecule", "resid-block"],
                   default="type-molecule",
                   help=(
                       "Chain grouping strategy.\n"
                       "  molnum        – mol column in dump (best if available)\n"
                       "  type-molecule – split polymer types into equal chains,\n"
                       "                  keep receptor types as single chains [DEFAULT]\n"
                       "  atomtype      – one group per atom type\n"
                       "  resid-block   – requires --residues-per-chain"
                   ))
    p.add_argument("--residues-per-chain", type=int, default=None,
                   help="Residues per chain (only for --chain-by resid-block)")
    p.add_argument("--polymer-chain-len", type=int, default=POLYMER_CHAIN_LEN,
                   help=f"Atoms per polymer chain for type-molecule mode (default: {POLYMER_CHAIN_LEN})")
    p.add_argument("--polymer-types", nargs="+", default=sorted(POLYMER_TYPES),
                   help="Atom type strings that are polymer beads (default: 1 2)")
    p.add_argument("--receptor-types", nargs="+", default=sorted(RECEPTOR_TYPES),
                   help="Atom type strings that are receptor beads (default: 3 4)")

    p.add_argument("--time-unit", choices=["ps", "ns", "us"], default="ps")
    p.add_argument("--movie-plane", choices=["xy", "xz", "yz"], default="xz",
                   help="Projection plane for HTML movie (default: xz — best for slab)")
    p.add_argument("--make-html-movie", action="store_true")
    p.add_argument("--movie-out", default=None)
    p.add_argument("--movie-stride", type=int, default=1)
    p.add_argument("--stride", type=int, default=1,
                   help="Analyse every Nth frame (default: 1)")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


# ═════════════════════════════════════════════════════════════════════════════
#  Chain grouping
# ═════════════════════════════════════════════════════════════════════════════

def get_chain_groups(u, chain_by, residues_per_chain=None,
                     polymer_chain_len=POLYMER_CHAIN_LEN,
                     polymer_types=None, receptor_types=None):
    """
    Returns a list of AtomGroups, one per chain.
    Receptor chains are always appended at the end so their 1-based indices
    equal (n_polymer_chains + 1), (n_polymer_chains + 2), ...
    """
    if polymer_types  is None: polymer_types  = sorted(POLYMER_TYPES)
    if receptor_types is None: receptor_types = sorted(RECEPTOR_TYPES)

    atoms = u.atoms
    types = np.asarray(atoms.types)

    # ── molnum ────────────────────────────────────────────────────────────────
    if chain_by == "molnum":
        if not hasattr(atoms, "molnums"):
            raise RuntimeError(
                "Topology has no molnums. Use --chain-by type-molecule instead."
            )
        molnums = np.asarray(atoms.molnums)
        uniq    = np.unique(molnums)
        groups  = [atoms[molnums == m] for m in uniq]
        return groups

    # ── atomtype (one group per type) ─────────────────────────────────────────
    if chain_by == "atomtype":
        uniq = np.unique(types)
        return [atoms[types == t] for t in uniq]

    # ── type-molecule (PRIMARY mode for this system) ──────────────────────────
    if chain_by == "type-molecule":
        groups = []

        # --- polymer chains: split each polymer type into equal blocks --------
        for pt in polymer_types:
            mask    = types == pt
            indices = np.where(mask)[0]
            n       = len(indices)
            if n == 0:
                continue
            if n % polymer_chain_len != 0:
                raise RuntimeError(
                    f"Atom type '{pt}' has {n} atoms, not divisible by "
                    f"polymer_chain_len={polymer_chain_len}. "
                    f"Adjust --polymer-chain-len."
                )
            for start in range(0, n, polymer_chain_len):
                groups.append(atoms[indices[start:start + polymer_chain_len]])

        # --- receptor chains: one group per receptor type ---------------------
        for rt in receptor_types:
            mask = types == rt
            if mask.any():
                groups.append(atoms[mask])

        return groups

    # ── resid-block ───────────────────────────────────────────────────────────
    if chain_by == "resid-block":
        if not residues_per_chain or residues_per_chain < 1:
            raise RuntimeError("--residues-per-chain required for resid-block")
        residues = u.residues
        n_res    = len(residues)
        if n_res % residues_per_chain != 0:
            raise RuntimeError(
                f"Total residues ({n_res}) not divisible by "
                f"residues_per_chain ({residues_per_chain}). "
                f"Try --chain-by type-molecule instead."
            )
        return [residues[i:i + residues_per_chain].atoms
                for i in range(0, n_res, residues_per_chain)]

    raise RuntimeError(f"Unknown --chain-by value: {chain_by}")


# ═════════════════════════════════════════════════════════════════════════════
#  Graph / geometry utilities
# ═════════════════════════════════════════════════════════════════════════════

def connected_components(adjacency: list[list[int]]) -> list[list[int]]:
    n    = len(adjacency)
    seen = np.zeros(n, dtype=bool)
    components = []
    for start in range(n):
        if seen[start]:
            continue
        queue, seen[start], comp = deque([start]), True, []
        while queue:
            node = queue.popleft()
            comp.append(node)
            for nbr in adjacency[node]:
                if not seen[nbr]:
                    seen[nbr] = True
                    queue.append(nbr)
        components.append(comp)
    return components


def minimum_image_vector(r1: np.ndarray, r2: np.ndarray,
                         box: np.ndarray) -> np.ndarray:
    dr      = r2 - r1
    lengths = box[:3]
    dr     -= lengths * np.round(dr / lengths)
    return dr


def chain_masses_and_coms(chain_groups):
    masses = np.zeros(len(chain_groups))
    coms   = np.zeros((len(chain_groups), 3))
    for i, ag in enumerate(chain_groups):
        m = None
        if hasattr(ag, "masses") and ag.masses is not None:
            m = np.asarray(ag.masses, dtype=float)
        total = float(m.sum()) if m is not None else 0.0
        if total > 0.0:
            masses[i] = total
            coms[i]   = ag.center_of_mass()
        else:
            masses[i] = float(len(ag))
            coms[i]   = ag.positions.mean(axis=0)
    return masses, coms


def unwrap_cluster_coms(cluster_indices, coms, box, adjacency):
    if len(cluster_indices) == 1:
        return np.array([coms[cluster_indices[0]]])
    cluster_set = set(cluster_indices)
    root        = cluster_indices[0]
    placed      = {root: coms[root].copy()}
    queue       = deque([root])
    while queue:
        i = queue.popleft()
        for j in adjacency[i]:
            if j not in cluster_set or j in placed:
                continue
            placed[j] = placed[i] + minimum_image_vector(coms[i], coms[j], box)
            queue.append(j)
    return np.array([placed[idx] for idx in cluster_indices])


def radius_of_gyration(points: np.ndarray, weights: np.ndarray) -> float:
    wsum = weights.sum()
    if wsum <= 0:
        return 0.0
    center = np.sum(points * weights[:, None], axis=0) / wsum
    rg2    = np.sum(weights * np.sum((points - center) ** 2, axis=1)) / wsum
    return float(np.sqrt(max(rg2, 0.0)))


def build_adjacency(selected_chain_groups, box, cutoff_angstrom,
                    min_pairs, backend):
    n         = len(selected_chain_groups)
    adjacency = [[] for _ in range(n)]
    for i in range(n - 1):
        pos_i = selected_chain_groups[i].positions
        if not len(pos_i):
            continue
        for j in range(i + 1, n):
            pos_j = selected_chain_groups[j].positions
            if not len(pos_j):
                continue
            dmat   = distance_array(pos_i, pos_j, box=box, backend=backend)
            npairs = int(np.count_nonzero(dmat <= cutoff_angstrom))
            if npairs >= min_pairs:
                adjacency[i].append(j)
                adjacency[j].append(i)
    return adjacency


# ═════════════════════════════════════════════════════════════════════════════
#  Receptor-specific analysis
# ═════════════════════════════════════════════════════════════════════════════

def receptor_wall_distance(chain_groups, receptor_idx_0based):
    """
    Min z-distance from each receptor chain to the nearest wall boundary
    (WALL_ZLO or WALL_ZHI).  Returns {1based_idx: dist_Ang}.
    """
    result = {}
    for idx in receptor_idx_0based:
        z       = chain_groups[idx].positions[:, 2]
        d_zlo   = np.abs(z - WALL_ZLO).min()
        d_zhi   = np.abs(z - WALL_ZHI).min()
        result[idx + 1] = float(min(d_zlo, d_zhi))
    return result


def receptor_contacts(adjacency, receptor_idx_0based):
    """
    Number of non-receptor chains within cutoff of each receptor.
    Returns {1based_idx: count}.
    """
    receptor_set = set(receptor_idx_0based)
    result = {}
    for idx in receptor_idx_0based:
        result[idx + 1] = sum(1 for nbr in adjacency[idx]
                              if nbr not in receptor_set)
    return result


def receptor_in_same_cluster(components, receptor_idx_0based):
    """True if all receptor chains share a single cluster."""
    rs = set(receptor_idx_0based)
    return any(rs.issubset(set(comp)) for comp in components)


# ═════════════════════════════════════════════════════════════════════════════
#  I/O helpers
# ═════════════════════════════════════════════════════════════════════════════

def write_two_column_dat(path, header, x, y):
    with Path(path).open("w") as fh:
        fh.write(header + "\n")
        for xi, yi in zip(x, y):
            fh.write(f"{xi:.6f} {yi:.10f}\n")


def plot_timeseries(x, y, xlabel, ylabel, title, outfile):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x, y, lw=1.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def plot_receptor_wall_distance(times_ps, wall_dist_series,
                                receptor_ids_1based, out_prefix):
    fig, ax = plt.subplots(figsize=(7, 5))
    for rid in receptor_ids_1based:
        ax.plot(times_ps, [d[rid] for d in wall_dist_series],
                lw=1.8, label=f"Receptor type {rid}")
    ax.axhline(WALL_CUTOFF, ls="--", color="red", lw=1.2,
               label=f"LJ cutoff ({WALL_CUTOFF} Å)")
    ax.axhline(0.0, ls=":", color="grey", lw=0.8)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Min distance to wall boundary (Å)")
    ax.set_title(
        "Receptor–wall distance vs time\n"
        f"wall/lj126  zlo={WALL_ZLO} Å  zhi={WALL_ZHI} Å  "
        f"ε={WALL_EPSILON}  σ={WALL_SIGMA}"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{out_prefix}receptor_wall_distance.png", dpi=300)
    plt.close(fig)


def plot_receptor_contacts(times_ps, contacts_series,
                           receptor_ids_1based, out_prefix):
    fig, ax = plt.subplots(figsize=(7, 5))
    for rid in receptor_ids_1based:
        ax.plot(times_ps, [c[rid] for c in contacts_series],
                lw=1.8, label=f"Receptor type {rid}")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Number of contacting chains")
    ax.set_title("Receptor contact count vs time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{out_prefix}receptor_contacts.png", dpi=300)
    plt.close(fig)


def plot_receptor_same_cluster(times_ps, same_cluster_series, out_prefix):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(times_ps, same_cluster_series.astype(float),
            lw=1.5, color="purple")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["separate", "same cluster"])
    ax.set_xlabel("Time (ps)")
    ax.set_title("Are both receptors in the same cluster?")
    fig.tight_layout()
    fig.savefig(f"{out_prefix}receptor_same_cluster.png", dpi=300)
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
#  HTML movie helpers
# ═════════════════════════════════════════════════════════════════════════════

def unwrap_points(points: np.ndarray, box_lengths: np.ndarray) -> np.ndarray:
    if len(points) <= 1:
        return points.copy()
    out    = np.zeros_like(points)
    out[0] = points[0]
    for i in range(1, len(points)):
        out[i] = out[0] + minimum_image_vector(points[0], points[i], box_lengths)
    return out


def cluster_centers_from_chain_coms(components, coms, box_lengths):
    centers = np.zeros((len(components), 3))
    sizes   = np.zeros(len(components), dtype=int)
    for cid, comp in enumerate(components):
        pts           = coms[np.array(comp, dtype=int)]
        centers[cid]  = unwrap_points(pts, box_lengths).mean(axis=0)
        sizes[cid]    = len(comp)
    return centers, sizes


def recenter_points(points: np.ndarray, box_lengths: np.ndarray) -> np.ndarray:
    if not len(points):
        return points.copy()
    ref     = points.mean(axis=0)
    shifted = points - ref
    shifted -= box_lengths * np.round(shifted / box_lengths)
    return shifted


def choose_plane(coords: np.ndarray, plane: str):
    if plane == "xy":
        return coords[:, 0], coords[:, 1], "x (Å)", "y (Å)"
    if plane == "xz":
        return coords[:, 0], coords[:, 2], "x (Å)", "z (Å)"
    return coords[:, 1], coords[:, 2], "y (Å)", "z (Å)"


def marker_sizes_from_population(populations: np.ndarray) -> np.ndarray:
    return 22.0 + 18.0 * np.sqrt(populations.astype(float))


def color_list(n: int):
    base = (pc.qualitative.Plotly + pc.qualitative.D3
            + pc.qualitative.G10  + pc.qualitative.T10
            + pc.qualitative.Alphabet)
    return [base[i % len(base)] for i in range(n)]


def write_html_movie(movie_frames, plane: str, html_out: str):
    max_clusters = max(fd["nclusters"] for fd in movie_frames)
    palette      = color_list(max_clusters)

    centered_all          = np.concatenate([fd["centers"] for fd in movie_frames])
    x_all, y_all, xl, yl = choose_plane(centered_all, plane)
    pad_x  = 0.10 * (x_all.max() - x_all.min() + 1e-6)
    pad_y  = 0.10 * (y_all.max() - y_all.min() + 1e-6)
    x_range = [float(x_all.min() - pad_x), float(x_all.max() + pad_x)]
    y_range = [float(y_all.min() - pad_y), float(y_all.max() + pad_y)]

    def make_trace(fd):
        x, y, _, _ = choose_plane(fd["centers"], plane)
        colors = [palette[i % len(palette)] for i in range(fd["nclusters"])]
        hover  = [
            f"cluster={cid}<br>pop={pop}<br>frac={frac:.1f}%"
            for cid, (pop, frac) in enumerate(
                zip(fd["sizes"], fd["fractions"]), start=1)
        ]
        return go.Scatter(
            x=x, y=y, mode="markers", text=hover, hoverinfo="text",
            marker=dict(
                size=marker_sizes_from_population(fd["sizes"]),
                color=colors, opacity=0.82,
                line=dict(width=1.0, color="black"),
            ),
        )

    plot_frames = [
        go.Frame(
            name=str(i),
            data=[make_trace(fd)],
            layout=go.Layout(
                title=(f"Cluster cartoon | time={fd['time_label']} | "
                       f"clusters={fd['nclusters']} | largest={fd['largest']}")
            ),
        )
        for i, fd in enumerate(movie_frames)
    ]

    initial = movie_frames[0]
    fig = go.Figure(
        data=[make_trace(initial)],
        layout=go.Layout(
            title=(f"Cluster cartoon | time={initial['time_label']} | "
                   f"clusters={initial['nclusters']} | largest={initial['largest']}"),
            xaxis=dict(title=xl, range=x_range, visible=False),
            yaxis=dict(title=yl, range=y_range, visible=False,
                       scaleanchor="x", scaleratio=1),
            template="plotly_white", showlegend=False,
            updatemenus=[dict(
                type="buttons", showactive=False,
                x=0.02, y=1.12, direction="left",
                buttons=[
                    dict(label="Play", method="animate",
                         args=[None, {"frame": {"duration": 100, "redraw": True},
                                      "transition": {"duration": 0},
                                      "fromcurrent": True}]),
                    dict(label="Pause", method="animate",
                         args=[[None], {"frame": {"duration": 0, "redraw": False},
                                        "mode": "immediate",
                                        "transition": {"duration": 0}}]),
                ],
            )],
            sliders=[dict(
                active=0, pad={"t": 35},
                steps=[
                    dict(method="animate",
                         label=f"{fd['time_value']:.2f}",
                         args=[[str(i)],
                               {"frame": {"duration": 0, "redraw": True},
                                "mode": "immediate",
                                "transition": {"duration": 0}}])
                    for i, fd in enumerate(movie_frames)
                ],
            )],
        ),
        frames=plot_frames,
    )
    fig.write_html(html_out, include_plotlyjs=True, full_html=True)


# ═════════════════════════════════════════════════════════════════════════════
#  main
# ═════════════════════════════════════════════════════════════════════════════

def main():
    args = parse_args()

    if args.cutoff_nm <= 0:
        raise ValueError("--cutoff-nm must be positive.")
    if args.min_pairs < 1:
        raise ValueError("--min-pairs must be >= 1.")
    if args.movie_stride < 1:
        raise ValueError("--movie-stride must be >= 1.")
    if args.stride < 1:
        raise ValueError("--stride must be >= 1.")

    cutoff_angstrom = args.cutoff_nm * 10.0
    tscale          = {"ps": 1.0, "ns": 1e-3, "us": 1e-6}[args.time_unit]
    out_prefix      = args.out_prefix          # default "receptor_"

    polymer_types  = [str(t) for t in args.polymer_types]
    receptor_types = [str(t) for t in args.receptor_types]

    # ── load universe ─────────────────────────────────────────────────────────
    if args.verbose:
        print(f"Loading topology : {args.topology}")
        print(f"Loading trajectory: {args.trajectory}")

    u = mda.Universe(
        args.topology,
        args.trajectory,
        topology_format="DATA",
        format="LAMMPSDUMP",
    )

    if args.verbose:
        types_u, counts_u = np.unique(u.atoms.types, return_counts=True)
        print("Atom type counts:", dict(zip(types_u, counts_u)))
        print(f"Total atoms: {len(u.atoms)}   residues: {len(u.residues)}")

    # ── build chain groups ────────────────────────────────────────────────────
    chain_groups = get_chain_groups(
        u,
        chain_by           = args.chain_by,
        residues_per_chain = args.residues_per_chain,
        polymer_chain_len  = args.polymer_chain_len,
        polymer_types      = polymer_types,
        receptor_types     = receptor_types,
    )
    n_chains = len(chain_groups)

    if n_chains < 2:
        raise RuntimeError(f"Only {n_chains} chain(s) found — need at least 2.")

    # ── identify receptor chain indices (0-based) ─────────────────────────────
    # In type-molecule mode receptors are appended last, one per receptor type.
    # Determine by checking which groups consist only of receptor-type atoms.
    all_types_set = set(receptor_types)
    receptor_idx_0based = []
    for i, ag in enumerate(chain_groups):
        grp_types = set(np.unique(np.asarray(ag.types)))
        if grp_types.issubset(all_types_set):
            receptor_idx_0based.append(i)

    if not receptor_idx_0based:
        raise RuntimeError(
            f"No receptor chains found for types {receptor_types}. "
            "Check --receptor-types."
        )

    receptor_ids_1based = [i + 1 for i in receptor_idx_0based]

    # ── atom selection for distance calculation ───────────────────────────────
    if args.atom_names:
        sel = " or ".join(f"name {n}" for n in args.atom_names)
        selected_chain_groups = [cg.select_atoms(sel) for cg in chain_groups]
        empty = [i for i, ag in enumerate(selected_chain_groups) if not len(ag)]
        if empty:
            raise RuntimeError(
                f"Chains {empty[:10]} have no atoms matching '{sel}'."
            )
    else:
        # use all atoms in each chain
        selected_chain_groups = chain_groups

    if args.verbose:
        print(f"\nChain-by mode    : {args.chain_by}")
        print(f"Total chains     : {n_chains}")
        print(f"Receptor chains  : indices (1-based) {receptor_ids_1based}")
        print(f"Cutoff           : {args.cutoff_nm:.3f} nm  ({cutoff_angstrom:.2f} Å)")
        print(f"Min pairs        : {args.min_pairs}")
        print(f"Wall zlo/zhi     : {WALL_ZLO}/{WALL_ZHI} Å")
        print(f"Output prefix    : {out_prefix}\n")

    # ── accumulators ──────────────────────────────────────────────────────────
    times_ps                 = []
    nclusters_series         = []
    largest_rg_series        = []
    largest_fraction_series  = []
    mean_cluster_size_series = []
    wall_dist_series         = []
    contacts_series          = []
    same_cluster_series      = []
    movie_frames             = []

    movie_out = (args.movie_out or f"{out_prefix}cluster_movie.html")

    cluster_dist_file    = Path(f"{out_prefix}cluster_size_distribution.dat")
    receptor_detail_file = Path(f"{out_prefix}receptor_details.dat")

    traj     = u.trajectory[::args.stride]
    iterator = (tqdm(traj, desc="Analysing frames", unit="frame")
                if args.verbose else traj)

    with (cluster_dist_file.open("w") as dist_fh,
          receptor_detail_file.open("w") as rec_fh):

        dist_fh.write("# time_ps cluster_size count\n")

        # receptor detail header
        rec_fh.write(
            "# time_ps  "
            + "  ".join(f"contacts_rectype{r}" for r in receptor_types)
            + "  "
            + "  ".join(f"wall_dist_rectype{r}_Ang" for r in receptor_types)
            + "  same_cluster\n"
        )

        for ts in iterator:
            time_ps = float(ts.time)
            if time_ps < args.tmin_ps:
                continue
            if args.tmax_ps is not None and time_ps > args.tmax_ps:
                break

            adjacency  = build_adjacency(
                selected_chain_groups, ts.dimensions,
                cutoff_angstrom, args.min_pairs, args.backend,
            )
            components = connected_components(adjacency)
            sizes      = np.array([len(c) for c in components], dtype=int)

            nclusters         = len(components)
            mean_cluster_size = float(sizes.mean())
            largest_idx       = int(np.argmax(sizes))
            largest_cluster   = components[largest_idx]
            largest_fraction  = 100.0 * len(largest_cluster) / n_chains

            masses, coms = chain_masses_and_coms(chain_groups)
            lc_coords    = unwrap_cluster_coms(
                largest_cluster, coms, ts.dimensions, adjacency)
            lc_masses    = masses[np.array(largest_cluster, dtype=int)]
            largest_rg   = radius_of_gyration(lc_coords, lc_masses) / 10.0

            # receptor metrics
            wall_dists = receptor_wall_distance(chain_groups, receptor_idx_0based)
            contacts   = receptor_contacts(adjacency, receptor_idx_0based)
            same_clust = receptor_in_same_cluster(components, receptor_idx_0based)

            # write receptor detail line
            rec_fh.write(
                f"{time_ps:.6f}  "
                + "  ".join(str(contacts[i + 1]) for i in receptor_idx_0based)
                + "  "
                + "  ".join(f"{wall_dists[i + 1]:.6f}" for i in receptor_idx_0based)
                + f"  {int(same_clust)}\n"
            )

            # cluster size distribution
            size_counts = Counter(sizes.tolist())
            for sz in sorted(size_counts):
                dist_fh.write(f"{time_ps:.6f} {sz:d} {size_counts[sz]:d}\n")

            times_ps.append(time_ps)
            nclusters_series.append(nclusters)
            largest_rg_series.append(largest_rg)
            largest_fraction_series.append(largest_fraction)
            mean_cluster_size_series.append(mean_cluster_size)
            wall_dist_series.append(wall_dists)
            contacts_series.append(contacts)
            same_cluster_series.append(same_clust)

            if args.make_html_movie and (len(times_ps) - 1) % args.movie_stride == 0:
                box_lengths       = np.asarray(ts.dimensions[:3])
                centers, csizes   = cluster_centers_from_chain_coms(
                    components, coms, box_lengths)
                centers   = recenter_points(centers, box_lengths)
                fractions = 100.0 * csizes / float(n_chains)
                dt        = time_ps * tscale
                movie_frames.append({
                    "time_value": dt,
                    "time_label": f"{dt:.3f} {args.time_unit}",
                    "centers":    centers.copy(),
                    "sizes":      csizes.copy(),
                    "fractions":  fractions.copy(),
                    "nclusters":  len(components),
                    "largest":    int(csizes.max()) if len(csizes) else 0,
                })

            if args.verbose and hasattr(iterator, "set_postfix"):
                iterator.set_postfix(
                    clusters  = nclusters,
                    largest   = f"{largest_fraction:.1f}%",
                    rg        = f"{largest_rg:.2f}nm",
                )

    # ── finalise arrays ───────────────────────────────────────────────────────
    times_ps                  = np.asarray(times_ps)
    nclusters_series          = np.asarray(nclusters_series)
    largest_rg_series         = np.asarray(largest_rg_series)
    largest_fraction_series   = np.asarray(largest_fraction_series)
    mean_cluster_size_series  = np.asarray(mean_cluster_size_series)
    same_cluster_series       = np.asarray(same_cluster_series, dtype=bool)

    if not times_ps.size:
        raise RuntimeError("No frames selected. Check --tmin-ps / --tmax-ps.")

    # ── write .dat files ──────────────────────────────────────────────────────
    write_two_column_dat(f"{out_prefix}nclusters.dat",
                         "# time_ps number_of_clusters",
                         times_ps, nclusters_series)
    write_two_column_dat(f"{out_prefix}largest_cluster_rg.dat",
                         "# time_ps largest_cluster_rg_nm",
                         times_ps, largest_rg_series)
    write_two_column_dat(f"{out_prefix}largest_cluster_fraction.dat",
                         "# time_ps largest_cluster_fraction_percent",
                         times_ps, largest_fraction_series)
    write_two_column_dat(f"{out_prefix}mean_cluster_size.dat",
                         "# time_ps mean_cluster_size",
                         times_ps, mean_cluster_size_series)

    # ── standard plots ────────────────────────────────────────────────────────
    plot_timeseries(times_ps, nclusters_series,
                    "Time (ps)", "Number of clusters",
                    "Number of clusters vs time",
                    f"{out_prefix}nclusters.png")
    plot_timeseries(times_ps, largest_rg_series,
                    "Time (ps)", "Largest cluster Rg (nm)",
                    "Largest cluster Rg vs time",
                    f"{out_prefix}largest_cluster_rg.png")
    plot_timeseries(times_ps, largest_fraction_series,
                    "Time (ps)", "Chains in largest cluster (%)",
                    "Largest cluster fraction vs time",
                    f"{out_prefix}largest_cluster_fraction.png")
    plot_timeseries(times_ps, mean_cluster_size_series,
                    "Time (ps)", "Mean cluster size",
                    "Mean cluster size vs time",
                    f"{out_prefix}mean_cluster_size.png")

    # ── receptor plots ────────────────────────────────────────────────────────
    plot_receptor_wall_distance(times_ps, wall_dist_series,
                                receptor_ids_1based, out_prefix)
    plot_receptor_contacts(times_ps, contacts_series,
                           receptor_ids_1based, out_prefix)
    plot_receptor_same_cluster(times_ps, same_cluster_series, out_prefix)

    # ── HTML movie ────────────────────────────────────────────────────────────
    if args.make_html_movie:
        if not movie_frames:
            raise RuntimeError("No movie frames collected.")
        write_html_movie(movie_frames, args.movie_plane, movie_out)

    # ── summary ───────────────────────────────────────────────────────────────
    if args.verbose:
        files = [
            f"{out_prefix}nclusters.dat",
            f"{out_prefix}nclusters.png",
            f"{out_prefix}largest_cluster_rg.dat",
            f"{out_prefix}largest_cluster_rg.png",
            f"{out_prefix}largest_cluster_fraction.dat",
            f"{out_prefix}largest_cluster_fraction.png",
            f"{out_prefix}mean_cluster_size.dat",
            f"{out_prefix}mean_cluster_size.png",
            f"{out_prefix}cluster_size_distribution.dat",
            f"{out_prefix}receptor_details.dat",
            f"{out_prefix}receptor_wall_distance.png",
            f"{out_prefix}receptor_contacts.png",
            f"{out_prefix}receptor_same_cluster.png",
        ]
        if args.make_html_movie:
            files.append(movie_out)
        print("\nWrote:")
        for f in files:
            print(f"  {f}")


if __name__ == "__main__":
    main()
Done

