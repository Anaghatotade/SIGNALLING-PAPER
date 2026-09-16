#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIM_BASE="$SCRIPT_DIR/../simulations/cutoff2.5_standard_LJ"

root_dirs=(
"$SIM_BASE/eps0.50_co_phase_separated"
"$SIM_BASE/eps0.25_distinct_condensates"
"$SIM_BASE/eps0.10_distinct_condensates"
"$SIM_BASE/eps0.40_touching_condensates"
)

SCRIPTS_DIR="$SCRIPT_DIR/scripts"
receptor_source_dir="$SIM_BASE/eps0.50_co_phase_separated/N_50"

for root in "${root_dirs[@]}"; do
    for dir in "$root"/N_*; do
        if [ -d "$dir" ]; then
            cp "$SCRIPTS_DIR/cluster_lammps.py" "$dir/"
            cp "$SCRIPTS_DIR/a_cluster_lammps.py" "$dir/"
            cp "$SCRIPTS_DIR/b_cluster_lammps.py" "$dir/"
            cp "$SCRIPTS_DIR/c_cluster_lammps.py" "$dir/"
            cp "$SCRIPTS_DIR/d_cluster_lammps.py" "$dir/"
            # receptor_cluster_lammps.py differs slightly between a few run
            # folders in the original data, so it is not deduplicated; this
            # copies the eps=0.50 version as the canonical template.
            cp "$receptor_source_dir/receptor_cluster_lammps.py" "$dir/"
            echo "Copied -> $dir"
        fi
    done
done
