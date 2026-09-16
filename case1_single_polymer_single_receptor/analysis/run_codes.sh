#!/bin/bash
# Runs a_cluster_lammps.py / b_cluster_lammps.py inside every N_* run folder.
# Requires that the (large, non-version-controlled) lammps.data and
# traj.lammpstrj files already exist in each run folder -- regenerate them
# by running the LAMMPS input scripts (see simulations/<condition>/N_*/lammps_sim.in)
# and converting the binary dump with `binary2txt file_trajectory.bin`.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

root_dirs=(
    "$SCRIPT_DIR/../simulations/LJ_eps0.5_phase_separation"
    "$SCRIPT_DIR/../simulations/WCA_no_phase_separation"
)

for root in "${root_dirs[@]}"; do
    echo "Processing $root"

    for dir in "$root"/N_*; do
        [ -d "$dir" ] || continue

        echo "=================================================="
        echo "Entering $dir"
        cd "$dir" || continue

        echo "Running a_cluster_lammps.py..."
        python a_cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 1 \
            --chain-by resid \
            --min-pairs 1 \
            --verbose \
            > a_cluster.log 2>&1

        echo "a_cluster_lammps.py DONE"

        echo "Running b_cluster_lammps.py..."
        python b_cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 2 \
            --chain-by resid \
            --min-pairs 1 \
            --verbose \
            > b_cluster.log 2>&1

        echo "b_cluster_lammps.py DONE"

        cd - >/dev/null
    done
done

echo "All jobs completed."
