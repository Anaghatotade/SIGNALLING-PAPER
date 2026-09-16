#!/bin/bash
# Deploys the canonical cluster-analysis scripts (analysis/scripts/) into every
# receptor-count run folder (N_*) under both simulation conditions, so that
# a_cluster_lammps.py / b_cluster_lammps.py can be run there against the
# (locally regenerated, not version-controlled) lammps.data / traj.lammpstrj
# trajectory files of that run.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/scripts"

root_dirs=(
    "$SCRIPT_DIR/../simulations/LJ_eps0.5_phase_separation"
    "$SCRIPT_DIR/../simulations/WCA_no_phase_separation"
)

for root in "${root_dirs[@]}"; do
    for dir in "$root"/N_*; do
        if [ -d "$dir" ]; then
            cp "$SOURCE_DIR/a_cluster_lammps.py" "$dir/"
            cp "$SOURCE_DIR/b_cluster_lammps.py" "$dir/"

            echo "Copied to: $dir"
        fi
    done
done

echo "All files copied successfully."
