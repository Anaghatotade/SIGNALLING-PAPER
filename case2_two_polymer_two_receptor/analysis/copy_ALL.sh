#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIM_BASE="$SCRIPT_DIR/../simulations/cutoff2.5_standard_LJ"

# Source directory containing the notebooks
# NOTE: each N_* folder already has its own executed copy of these notebooks
# with run-specific saved output. Only use this script to seed a *fresh*,
# not-yet-run N_* folder with the template notebooks -- rerunning it against
# folders that already contain completed analysis will overwrite those results.
src_dir="$SCRIPT_DIR/scripts"

root_dirs=(
"$SIM_BASE/eps0.50_co_phase_separated"
"$SIM_BASE/eps0.25_distinct_condensates"
"$SIM_BASE/eps0.10_distinct_condensates"
"$SIM_BASE/eps0.40_touching_condensates"
)

# Notebook names
files=(
"largest_conn_cluster_receptorC.ipynb"
"largest_conn_cluster_receptorD.ipynb"
)

for root in "${root_dirs[@]}"; do
    echo "Processing: $root"

    for dir in "$root"/N_*; do
        [ -d "$dir" ] || continue

        echo "  -> Copying to $(basename "$dir")"

        for file in "${files[@]}"; do
            cp -f "$src_dir/$file" "$dir/"
        done
    done
done

echo "Done."
