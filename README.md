# Phase Separation of Intracellular Proteins Influences Extracellular Signaling at the Synapse

This repository contains the coarse-grained Langevin dynamics simulations, LAMMPS
input files, and analysis code used in:

> V. Saxena\*, A. A. Totade\*, V. Pandey, G. Chauhan, *"Phase Separation of
> Intracellular Proteins Influences Extracellular Signaling at the Synapse"*
> (\*equal contribution), Department of Chemical Engineering, Indian Institute
> of Technology Indore. Corresponding author: gaurav@iiti.ac.in

## Scientific objective

Neuronal signaling depends on receptor clustering at the plasma membrane, a
process influenced by intracellular scaffold proteins that can themselves
undergo liquid–liquid phase separation (LLPS) into membraneless condensates.
This project asks two questions with a minimal coarse-grained model:

1. **How does receptor clustering at the membrane depend on the phase
   behaviour (dispersed vs. condensed) of a bulk scaffold polymer that binds
   the receptors?**
2. **When two competing scaffold polymers are present, how does their
   relative spatial organization (separate condensates, touching
   condensates, or a single co-phase-separated condensate) modulate
   receptor clustering?**

The model couples a bead-spring (Kremer–Grest) scaffold polymer in a bulk
compartment to monomeric receptor beads confined to an adjacent membrane
compartment, connected through attractive polymer–receptor interactions
across a repulsive wall. See the paper (Sections I–III) for the full
scientific motivation, results, and discussion.

## Repository structure

The two simulation systems described in the paper are kept in separate,
self-contained folders:

```
SIGNALLING-PAPER/
├── case1_single_polymer_single_receptor/   # Section II.A/II.B & Fig. 2–3
└── case2_two_polymer_two_receptor/         # Section II.C & Fig. 4–5
```

Each case folder follows the same internal layout:

```
caseN_.../
├── README.md                  # case-specific methodology and file guide
├── initial_configuration/     # equilibration stage: builds the starting polymer configuration
├── simulations/                # production LAMMPS runs, organized by physical condition, then by receptor count (N_<R>)
├── analysis/                   # analysis scripts, notebooks, and results
└── exploratory/  (case 2 only) # trial runs not part of the paper's main published cases
```

See `case1_single_polymer_single_receptor/README.md` and
`case2_two_polymer_two_receptor/README.md` for the full breakdown of
conditions, receptor counts, and scripts in each case.

## The two simulation systems

### Case 1 — single scaffold polymer + single receptor species (Section II.A–B)

One scaffold polymer (Kremer–Grest bead-spring chain, N = 25 monomers) in a
bulk compartment, and a single receptor species (monomeric beads) confined to
an adjacent membrane compartment. Two scaffold-polymer interaction models are
compared at receptor counts R = 0, 1, 20, 50, 100, 125, 150, 200, 250, 300:

| Condition | Polymer (P–P) potential | Behaviour |
|---|---|---|
| `WCA_no_phase_separation` | purely repulsive WCA (cutoff 2¹ᐟ⁶σ) | scaffold does **not** phase separate |
| `LJ_eps0.5_phase_separation` | attractive LJ, ε = 0.5 kBT (cutoff 2.5σ) | scaffold **condenses** into a dense phase |

This corresponds to Fig. 2 and Fig. 3 of the paper: condensation of the
scaffold polymer lowers the critical receptor concentration for clustering by
~70% and produces larger receptor clusters at saturation.

### Case 2 — two scaffold polymers + two receptor species (Section II.C, Results III.B)

Two scaffold polymers (P1, P2), each with strong self-attraction
(P1–P1, P2–P2, ε = 0.5), and two receptor species (R1, R2) that each bind
strongly to only one of the two polymers (R1↔P1, R2↔P2). The inter-polymer
interaction strength ε(P1–P2) is varied to change how the two condensed
phases are spatially organized, at receptor counts R1 = R2 = 50, 100, 204:

| Condition (`simulations/cutoff2.5_standard_LJ/`) | ε(P1–P2) | Topology |
|---|---|---|
| `eps0.10_distinct_condensates` | 0.10 | two fully phase-separated condensates |
| `eps0.25_distinct_condensates` | 0.25 | two fully phase-separated condensates |
| `eps0.40_touching_condensates` | 0.40 | condensates just touching |
| `eps0.50_co_phase_separated`   | 0.50 | single co-phase-separated condensate |

This corresponds to Fig. 4 and Fig. 5 of the paper: co-phase separation of
the two polymers (ε = 0.5) *suppresses* receptor clustering relative to
distinct condensates, indicating that competing scaffold–scaffold
interactions can antagonize receptor clustering.

`case2_two_polymer_two_receptor/exploratory/` holds two additional sets of
runs that are **not** among the four published P1–P2 cases above: a
WCA-based (non-attractive) inter-polymer interaction control, and a set of
"core–shell" topology trial runs. Both are preserved for completeness but
are explicitly skipped by the paper's own analysis scripts (see the comment
in `analysis/rdf_e_basecase.py`).

## LAMMPS methodology (summary; see paper Sec. II for full detail)

- Implicit-solvent **Langevin dynamics**, integrated with the velocity-Verlet
  algorithm, damping parameter 1/τ, timestep 0.005τ, in reduced
  Lennard-Jones units.
- Scaffold polymers: Kremer–Grest bead-spring chains, N = 25 monomers per
  chain, bonded via a FENE potential (K = 10, R₀ = 1.0) and non-bonded via a
  Lennard-Jones potential (Eq. 1 in the paper), or a WCA potential (cutoff
  2¹ᐟ⁶σ) for the non-phase-separating control.
  See Table I of the paper for the full non-bonded parameter set (ε, σ,
  cutoff) between every pair of particle types (P1, P2, R1, R2).
- The simulation box is split by a repulsive WCA wall along z into a lower
  (scaffold-polymer) and upper (receptor) compartment, confining receptors
  to the membrane-proximal region and preventing polymers from crossing
  into it.
- Simulations were run with **LAMMPS** (23 Jun 2022 build) in the NVT
  ensemble for 2×10⁷ production steps; the last 200 equilibrated trajectory
  frames were used to compute cluster statistics. Trajectories were
  visualized with **OVITO**.

## Analysis workflow

1. **Initial configuration** (`initial_configuration/`): builds/equilibrates
   the starting scaffold-polymer configuration (`Input_file.ipynb` +
   `lammps_sim.in`), producing the restart file that the production runs
   read from.
2. **Production runs** (`simulations/<condition>/N_<R>/`): each folder holds
   the LAMMPS input script (`lammps_sim.in`) for that specific
   condition/receptor-count combination, plus its own equilibration log
   (`log.polymer_wall.txt`) and, once run, the (large, not version
   controlled) raw trajectory.
3. **Cluster analysis** (`analysis/`): standalone CLI tools
   (`a/b/[c/d]_cluster_lammps.py`, `cluster_lammps.py`,
   `receptor_cluster_lammps.py`) compute connected-component cluster
   statistics from a run's `lammps.data` + `traj.lammpstrj` files. The
   `copy_*.sh` helper scripts deploy these tools (and companion notebooks)
   into every `N_*` folder; the `run_*.sh` scripts then execute them.
4. **Aggregation** (`analysis/fraction_vs_N.py`,
   `analysis/largest_cluster_analysis.py`, `analysis/rdf*.py`,
   `analysis/collect_cluster_plots.py` [case 2 only]): walk every
   condition/receptor-count folder, compute the ensemble-averaged fraction
   of receptors/polymers in the largest cluster as a function of receptor
   count, and produce the summary plots analogous to Fig. 3 and Fig. 5 of
   the paper. Outputs land in `analysis/results/`.

See each case's own `README.md` for the exact command-line usage of every
script.

## Reproducing a simulation

```bash
cd case1_single_polymer_single_receptor/simulations/WCA_no_phase_separation/N_100
lmp -in lammps_sim.in          # requires a local LAMMPS build
```

Raw trajectories (`lammps.data`, `traj.lammpstrj`, restart checkpoints, and
`data*.xyz` OVITO-visualization snapshots) are regenerated by the run and are
**not** version-controlled (see `.gitignore`) because of their size. What is
version-controlled is everything needed to regenerate them: the LAMMPS input
scripts, the shared starting configuration, and the small numeric/plot
results already computed from past runs.


elsewhere; running it against this already-cleaned repository is a no-op
(all source paths will already have been moved).
