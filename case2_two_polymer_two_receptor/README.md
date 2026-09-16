# Case 2 — Two Scaffold Polymers + Two Receptor Species

Corresponds to Section II.C and Results III.B (Fig. 4, Fig. 5) of the paper.
This system asks: *when two scaffold polymers compete or co-condense, how
does their relative spatial organization modulate receptor clustering?*

## Model summary

- Two scaffold polymer species, P1 and P2, each a Kremer–Grest bead-spring
  chain (as in Case 1), each with strong self-attraction
  (P1–P1, P2–P2, ε = 0.5).
- Two receptor species, R1 and R2, confined to the membrane compartment.
  Each receptor binds strongly to only one polymer and is strongly repelled
  by the other: R1↔P1 attraction / R1↔P2 repulsion, and vice versa for R2.
- The inter-polymer interaction strength ε(P1–P2) is varied
  (0.10, 0.25, 0.40, 0.50) to change how the two polymers' condensed phases
  are organized relative to one another, at matched receptor counts
  R1 = R2 = 50, 100, 204.
- Box: 80 × 80 × 70.

## Directory layout

```
case2_two_polymer_two_receptor/
├── initial_configuration/          # (preliminary) equilibration stage; see note below
│   ├── Input_file.ipynb
│   ├── initial.xyz, lammps_sim.in, prod.restart1, log.polymer_wall.txt, submit_active.lsf
│
├── simulations/cutoff2.5_standard_LJ/     # the four published P1-P2 interaction cases (standard LJ cutoff = 2.5σ)
│   ├── eps0.10_distinct_condensates/       # ε(P1-P2)=0.10 -> two fully separated condensates
│   ├── eps0.25_distinct_condensates/       # ε(P1-P2)=0.25 -> two fully separated condensates
│   ├── eps0.40_touching_condensates/       # ε(P1-P2)=0.40 -> condensates just touching
│   ├── eps0.50_co_phase_separated/         # ε(P1-P2)=0.50 -> single co-phase-separated condensate
│   └── each condition/
│       ├── submit_N.lsf, submit_jup.lsf     # one shared HPC submission template for this condition
│       └── N_<R>/                            # one folder per (R1=R2=R) receptor count
│           ├── initial.xyz                   # full two-polymer starting configuration (identical across R and ε; see note below)
│           ├── lammps_sim.in                 # production input (reads initial.xyz)
│           ├── lammps_sim1.in                # continuation input (reads prod.restart1) for extending a run past one job's walltime
│           ├── A/B/C/D_*.dat, *.png           # per-species (polymerA, polymerB, receptorC, receptorD) cluster/RDF outputs
│           ├── largest_conn_cluster_polymerA.ipynb / polymerB.ipynb   # per-run notebooks (with saved results)
│           ├── receptor_cluster_lammps.py     # per-run copy of the CLI receptor-cluster tool (not fully identical across all runs; see note below)
│           └── plots.ipynb                    # per-run plotting notebook
│
├── exploratory/                     # NOT part of the paper's four published P1-P2 cases; kept for completeness
│   ├── wca_interpolymer_eps0.5_cutoff1.1225/   # P1-P2 interaction modeled as purely repulsive WCA instead of attractive LJ
│   └── core_shell_trials/try0 .. try5/          # preliminary "core-shell" topology trial runs
│
└── analysis/                        # case-wide analysis scripts & aggregated results
    ├── scripts/                      # canonical copies of a/b/c/d_cluster_lammps.py, cluster_lammps.py, largest_conn_cluster_receptorC/D.ipynb
    ├── copy_cluster.sh, run_abcluster.sh, run_cluster.sh   # deploy + run the CLI cluster-analysis tools across every N_*
    ├── copy_ALL.sh, run_ALL.sh                              # deploy + run the per-species notebooks across every N_*
    ├── copy_plots.sh, run_plots.sh                          # deploy + run plots.ipynb across every N_*
    ├── collect_cluster_plots.py                             # standalone MDAnalysis-based cluster-fraction-vs-frame plotter, collects results into analysis/results/all_largest_cluster_fraction_png/
    ├── fraction_vs_N.py, largest_cluster_analysis.py, rdf.py, rdf_e_basecase.py   # CLI scripts aggregating across the four eps(P1-P2) conditions -> Fig. 5-style plots
    └── results/                                             # all small numeric/plot outputs (csv, png, pdf)
```

## Important parameters (Table I of the paper)

| Particles | ε | σ | cutoff | Interaction |
|---|---|---|---|---|
| P1–P1 | 0.5 | 1 | 2.5σ | Attraction |
| P2–P2 | 0.5 | 1 | 2.5σ | Attraction |
| P1–P2 | 0.1 / 0.25 / 0.4 / 0.5 | 1 | 2.5σ | Weak attraction (varied — the four published cases) |
| R1–P1, R2–P2 | 50 | 1 | 4.0σ | Strong attraction |
| R1–P2, R2–P1 | 50 | 1 | 2¹ᐟ⁶σ | Strong repulsion |
| R1–R1, R1–R2, R2–R2 | 1.0 | 1 | 2¹ᐟ⁶σ | Repulsion |

## Notes on the source data (kept for transparency, not "fixed")

- **`initial.xyz` duplication.** The per-run `initial.xyz` files are
  byte-identical (30,000 atoms, polymer-only) across every receptor count
  (50/100/204) and every ε(P1–P2) condition, including the exploratory
  cases. They were left duplicated in place rather than deduplicated,
  because the checked-in `lammps_sim.in` for these runs contains no visible
  `create_atoms`/`delete_atoms` step that would explain how the receptor
  populations (R1, R2) differing between runs end up in the system; the
  full receptor-insertion step used in the original runs is therefore not
  fully reconstructable from the committed input scripts alone. Rather than
  guess at a missing preprocessing step, the original per-run files are
  preserved as-is.
- **`initial_configuration/`** here holds a preliminary equilibration
  (its `prod.restart1` differs from the one embedded in most of the
  production `N_*` folders), rather than being the direct one-to-one source
  of the production `initial.xyz`. It is kept for provenance but should not
  be assumed to be the literal precursor of every production run.
- **`receptor_cluster_lammps.py`** was *not* deduplicated to a single
  canonical script (unlike `a/b/c/d_cluster_lammps.py`, which were
  byte-identical everywhere) because a small number of run folders contain a
  slightly different version of it.
- **Exploratory folders.** `exploratory/wca_interpolymer_eps0.5_cutoff1.1225/`
  and `exploratory/core_shell_trials/` are additional simulation attempts
  present in the original repository. They are explicitly skipped by the
  paper's own `rdf_e_basecase.py` aggregation script (see the `# Will be
  skipped` comment) and are not part of the four main P1–P2 cases discussed
  in the paper's Results section. They are kept, clearly separated, for
  completeness rather than deleted.

## Analysis workflow

1. Prepare/equilibrate a starting configuration (see notes above).
2. Run a production simulation for a given ε(P1–P2) condition and receptor
   count:
   ```bash
   cd simulations/cutoff2.5_standard_LJ/eps0.50_co_phase_separated/N_100
   lmp -in lammps_sim.in
   ```
   Use `lammps_sim1.in` (`read_restart prod.restart1`) to continue a run
   from its own checkpoint if it needs to be split across multiple jobs.
3. Convert/prepare the resulting trajectory as `lammps.data` +
   `traj.lammpstrj` in the run folder (large, not version-controlled).
4. Run cluster analysis with the CLI tools in `analysis/scripts/`
   (`cluster_lammps.py` for combined polymer clusters, `a/c_cluster_lammps.py`
   for polymer A / receptor C individually, `b/d_cluster_lammps.py` for
   polymer B / receptor D, `receptor_cluster_lammps.py` for combined
   receptor clusters), deployed via `analysis/copy_cluster.sh` and executed
   with `analysis/run_abcluster.sh` / `run_cluster.sh`.
5. Aggregate across all four ε(P1–P2) conditions and receptor counts with
   `analysis/fraction_vs_N.py`, `analysis/largest_cluster_analysis.py`, and
   `analysis/rdf*.py` — reproducing the Fig. 5-style comparison of receptor
   clustering vs. inter-polymer interaction strength. Outputs land in
   `analysis/results/`.

## What each major script does

| Script | Purpose |
|---|---|
| `simulations/.../N_<R>/lammps_sim.in` | Production LAMMPS run for a given ε(P1–P2) and receptor count. |
| `simulations/.../N_<R>/lammps_sim1.in` | Continuation run from `prod.restart1` (same pair coefficients). |
| `analysis/scripts/cluster_lammps.py` | Connected-component clustering across both polymer species combined. |
| `analysis/scripts/a_cluster_lammps.py` / `c_cluster_lammps.py` | Clustering restricted to polymer A / receptor C (atom type 1 / type 3). |
| `analysis/scripts/b_cluster_lammps.py` / `d_cluster_lammps.py` | Clustering restricted to polymer B / receptor D (atom type 2 / type 4). |
| `receptor_cluster_lammps.py` (per-run) | Combined R1+R2 receptor clustering. |
| `analysis/collect_cluster_plots.py` | Standalone MDAnalysis-based largest-cluster-fraction-vs-frame calculator; collects one plot per run into `analysis/results/all_largest_cluster_fraction_png/`. |
| `analysis/fraction_vs_N.py` | Aggregates the largest-cluster fraction across every condition/receptor-count folder and plots vs. receptor count. |
| `analysis/largest_cluster_analysis.py` | OVITO-based cluster-size analysis across the four eps(P1-P2) conditions. |
| `analysis/rdf.py` / `rdf_e_basecase.py` | Radial distribution function comparison across conditions (the `_e_basecase` variant explicitly skips the exploratory WCA case). |
| `analysis/copy_cluster.sh` / `run_abcluster.sh` / `run_cluster.sh` | Deploy and execute the CLI cluster-analysis tools across every `N_*` folder. |
| `analysis/copy_ALL.sh` / `run_ALL.sh` | Deploy and execute the per-species result notebooks across every `N_*` folder. |
| `analysis/copy_plots.sh` / `run_plots.sh` | Deploy and execute `plots.ipynb` across every `N_*` folder. |
