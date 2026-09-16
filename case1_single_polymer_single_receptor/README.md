# Case 1 — Single Scaffold Polymer + Single Receptor Species

Corresponds to Sections II.A–II.B and Results III.A (Fig. 2, Fig. 3) of the
paper. This system asks: *does phase separation of a single bulk scaffold
polymer change the propensity of a membrane receptor to cluster?*

## Model summary

- One scaffold polymer species (P): a Kremer–Grest bead-spring chain,
  N = 25 monomers, FENE-bonded (K = 10, R₀ = 1.0).
- One receptor species (R): monomeric beads confined to a membrane
  compartment separated from the polymer's bulk compartment by a repulsive
  WCA wall.
- Two polymer self-interaction models are compared:
  - **WCA** (purely repulsive, cutoff 2¹ᐟ⁶σ) → polymer does not condense.
  - **LJ, ε = 0.5 kBT** (cutoff 2.5σ) → polymer condenses into a dense phase.
- Receptor count R (folder name `N_<R>`) swept over
  R = 0, 1, 20, 50, 100, 125, 150, 200, 250, 300 for the WCA condition, and
  R = 0, 1, 20, 50, 100, 150, 200, 250, 300 for the LJ ε=0.5 condition.
- Box: 60 × 60 × 70 (single polymer + single receptor per initial setup).

## Directory layout

```
case1_single_polymer_single_receptor/
├── initial_configuration/          # builds & equilibrates the starting polymer configuration
│   ├── Input_file.ipynb            # writes initial.xyz (fixed random seed 42)
│   ├── initial.xyz                 # starting coordinates (input to the equilibration run)
│   ├── lammps_sim.in               # equilibration LAMMPS script (minimize + Langevin dynamics)
│   ├── prod.restart1               # equilibrated restart checkpoint, read by every production run below
│   ├── log.polymer_wall.txt        # thermo log of the equilibration run
│   └── submit_active.lsf           # example HPC (LSF) batch submission script
│
├── simulations/
│   ├── WCA_no_phase_separation/    # polymer P-P = purely repulsive WCA (no condensation)
│   │   ├── submit_N.lsf, submit_jup.lsf   # one shared HPC submission template for this condition
│   │   ├── run_all_lammps.sh              # runs LAMMPS sequentially across all N_* folders locally
│   │   └── N_<R>/                          # one folder per receptor count R
│   │       ├── lammps_sim.in               # production input: reads prod.restart1, adds R receptors, runs 2e7 steps
│   │       ├── prod.restart1                # copy of the equilibrated starting configuration (see note below)
│   │       ├── log.polymer_wall.txt         # thermo log for this run
│   │       ├── a_cluster_lammps.py / b_cluster_lammps.py   # per-run copies of the cluster-analysis CLI tools
│   │       ├── largest_conn_cluster.ipynb, plots.ipynb, z_dense_phase_density_3_11.ipynb  # per-run analysis notebooks (with saved results)
│   │       └── Total_dense_chains.npy, dr_zero_positions.json  # small cached intermediate analysis arrays
│   │
│   └── LJ_eps0.5_phase_separation/ # polymer P-P = attractive LJ, ε=0.5 (condenses)
│       └── ... (same internal layout as above)
│
├── analysis/                        # case-wide analysis scripts & aggregated results
│   ├── scripts/                     # canonical copies of a_cluster_lammps.py, b_cluster_lammps.py
│   ├── copy_codes.sh, run_codes.sh          # deploy + run the CLI cluster-analysis tools across every N_*
│   ├── copy_plots.sh, run_plots.sh          # deploy + run plots.ipynb across every N_*
│   ├── Analysis.sh, run_all_notebooks.sh    # deploy + run the density-profile notebook across every N_*
│   ├── condensed_chain_fraction_vs_receptor_count.ipynb  # fraction of chains in the dense phase vs. R (uses Total_dense_chains.npy)
│   ├── z_dense_phase_density.ipynb          # z-density profile of the condensate
│   ├── radial_density_analysis.ipynb        # radial distribution / minimum-image-distance analysis
│   ├── largest_conn_cluster.ipynb           # standalone largest-cluster notebook (case-level copy)
│   ├── fraction_vs_N.py, largest_cluster_analysis.py, rdf_e_basecase.py   # CLI scripts aggregating across both conditions -> Fig. 3-style plots
│   └── results/                             # all small numeric/plot outputs (csv, png, pdf, npy caches)
│
└── NOTES_from_original_repo.txt     # short author notes recovered verbatim from the original repository
```

## Important parameters (Table I of the paper, subset relevant to this case)

| Particles | ε | σ | cutoff | Interaction |
|---|---|---|---|---|
| P–P (attractive case) | 0.5 | 1 | 2.5σ | Attraction (drives condensation) |
| P–P (WCA case) | 0.5 | 1 | 2¹ᐟ⁶σ | Purely repulsive (no condensation) |
| R–P | 50 | 1 | 4.0σ | Strong attraction (receptor binds scaffold) |
| R–R | 1.0 | 1 | 2¹ᐟ⁶σ | Repulsion |

From `NOTES_from_original_repo.txt`: linker length 2σ; walls placed at
z = ±22 and z = 25; in the eps_0.5 condition, receptors ("ligands") do not
interact with one another, only with the scaffold macromolecule.

## Note on `prod.restart1` duplication

Every `N_*` folder contains its own copy of `prod.restart1`, the restart
checkpoint written at the end of the equilibration stage in
`initial_configuration/`. In most run folders this file is byte-identical to
`initial_configuration/prod.restart1`; a few run folders in the original data
contain a restart checkpoint with a different hash. Because of this
inconsistency in the source data, these per-run copies were **not**
deduplicated during cleanup (to avoid silently altering which starting
configuration a given run actually used) — they are kept as-is for
provenance. If you need the canonical starting configuration, use
`initial_configuration/prod.restart1`.

## Analysis workflow

1. Equilibrate the scaffold polymer: run `initial_configuration/lammps_sim.in`
   (reads `initial.xyz`, writes `prod.restart1`).
2. Run a production simulation for a given condition/receptor count:
   ```bash
   cd simulations/WCA_no_phase_separation/N_100
   lmp -in lammps_sim.in
   ```
   This reads `prod.restart1`, adds R receptors via `create_atoms`, and runs
   2×10⁷ steps, appending to `log.polymer_wall.txt`.
3. Convert/prepare the resulting trajectory as `lammps.data` +
   `traj.lammpstrj` in the run folder (these raw outputs are large and are
   **not** version-controlled; see the repository-level `.gitignore`).
4. Run cluster analysis in that folder with the CLI tools in
   `analysis/scripts/` (deployed via `analysis/copy_codes.sh` then
   `analysis/run_codes.sh`), or open the per-run Jupyter notebooks directly.
5. Aggregate across all receptor counts and both conditions with
   `analysis/fraction_vs_N.py`, `analysis/largest_cluster_analysis.py`, and
   `analysis/rdf_e_basecase.py` — these reproduce the Fig. 3-style
   largest-cluster-fraction-vs-receptor-count comparison. Outputs are written
   to `analysis/results/`.

## What each major script does

| Script | Purpose |
|---|---|
| `initial_configuration/Input_file.ipynb` | Generates `initial.xyz`, the starting polymer-only configuration (fixed seed). |
| `simulations/<condition>/N_<R>/lammps_sim.in` | Production LAMMPS run: reads the equilibrated restart, inserts R receptors, integrates 2×10⁷ Langevin dynamics steps. |
| `analysis/scripts/a_cluster_lammps.py` / `b_cluster_lammps.py` | Connected-component cluster analysis CLI tools (operate on `lammps.data` + `traj.lammpstrj`; different atom-type/chain-grouping options select polymer vs. receptor clustering). |
| `analysis/condensed_chain_fraction_vs_receptor_count.ipynb` | Computes the fraction of scaffold chains residing in the dense phase as a function of receptor count. |
| `analysis/z_dense_phase_density.ipynb` | Computes the density profile of the condensate along z. |
| `analysis/radial_density_analysis.ipynb` | Radial density / minimum-image-distance analysis around the condensate. |
| `analysis/fraction_vs_N.py` | Aggregates the largest-cluster fraction across every `N_*` folder for both conditions and plots it vs. receptor count (Fig. 3 style). |
| `analysis/largest_cluster_analysis.py` | OVITO-based cluster-size analysis across all runs, with automatic interaction-cutoff detection from the folder path. |
| `analysis/rdf_e_basecase.py` | Radial distribution function comparison across conditions. |
| `analysis/copy_codes.sh` / `run_codes.sh` | Deploy and execute the CLI cluster-analysis tools across every `N_*` folder. |
| `analysis/copy_plots.sh` / `run_plots.sh` | Deploy and execute `plots.ipynb` across every `N_*` folder. |
| `analysis/Analysis.sh` / `run_all_notebooks.sh` | Deploy and execute the density-profile notebook across every `N_*` folder (`Analysis.sh` also submits each as an HPC batch job). |
