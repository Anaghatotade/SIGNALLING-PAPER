#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIM_BASE="$SCRIPT_DIR/../simulations/cutoff2.5_standard_LJ"

root_dirs=(
"$SIM_BASE/eps0.50_co_phase_separated"
"$SIM_BASE/eps0.25_distinct_condensates"
"$SIM_BASE/eps0.10_distinct_condensates"
"$SIM_BASE/eps0.40_touching_condensates"
)

for root in "${root_dirs[@]}"
do
    find "$root" -type f -name "cluster_lammps.py" | while read pyfile
    do
        workdir=$(dirname "$pyfile")

        echo "Running in: $workdir"

        (
            cd "$workdir" || exit

            python receptor_cluster_lammps.py \
                --lammps-data lammps.data \
                --trajectory traj.lammpstrj \
                --cutoff 2.5 \
                --atom-types 3 4 \
                --chain-by atom-block \
                --atoms-per-chain 200 \
                --min-pairs 1 \
                --verbose
                > receptor_cluster.log 2>&1
        )
    done
done
