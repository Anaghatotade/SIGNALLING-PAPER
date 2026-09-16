#!/bin/bash
set -uo pipefail   
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIM_BASE="$SCRIPT_DIR/../simulations/cutoff2.5_standard_LJ"

root_dirs=(
"$SIM_BASE/eps0.50_co_phase_separated"
"$SIM_BASE/eps0.25_distinct_condensates"
"$SIM_BASE/eps0.10_distinct_condensates"
"$SIM_BASE/eps0.40_touching_condensates"
)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

run_block() {
    local dir="$1"
    log "STARTING DIR: $dir"
    (
        set -e   
        cd "$dir"

        log "  Running cluster_lammps.py"
        python cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 1 2 \
            --chain-by atom-block \
            --atoms-per-chain 200 \
            --min-pairs 1 \
            --verbose \
            > polymer_cluster.log 2>&1
        log "  cluster_lammps.py DONE"

        log "  Running a_cluster_lammps.py (atom-type 1)"
        python c_cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 1 \
            --chain-by resid \
            --min-pairs 1 \
            --verbose \
            > a_cluster.log 2>&1
        log "  a_cluster_lammps.py DONE"

        log "  Running b_cluster_lammps.py"
        python b_cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 2 \
            --chain-by resid \
            --min-pairs 1 \
            --verbose \
            > b_cluster.log 2>&1
        log "  b_cluster_lammps.py DONE"

        log "  Running receptor_cluster_lammps.py"
        python receptor_cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 3 4 \
            --chain-by atom-block \
            --atoms-per-chain 200 \
            --min-pairs 1 \
            --verbose \
            > receptor_cluster.log 2>&1
        log "  receptor_cluster_lammps.py DONE"

        log "  Running c_cluster_lammps.py (atom-type 3)"
        python c_cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 3 \
            --chain-by resid \
            --min-pairs 1 \
            --verbose \
            > c_cluster.log 2>&1
        log "  c_cluster_lammps.py DONE"

        log "  Running d_cluster_lammps.py"
        python d_cluster_lammps.py \
            --lammps-data lammps.data \
            --trajectory traj.lammpstrj \
            --cutoff 2.5 \
            --atom-types 4 \
            --chain-by resid \
            --min-pairs 1 \
            --verbose \
            > d_cluster.log 2>&1
        log "  d_cluster_lammps.py DONE"
    )
    local status=$?

    if [ "$status" -eq 0 ]; then
        log "FINISHED SUCCESS: $dir"
    else
        log "FAILED (exit=$status): $dir"
    fi
    return "$status"
}


all_dirs=()
for root in "${root_dirs[@]}"; do
    if [ ! -d "$root" ]; then
        log "WARNING: root dir not found, skipping: $root"
        continue
    fi
    for dir in "$root"/N_*; do
        [ -d "$dir" ] || continue
        all_dirs+=("$dir")
    done
done
TOTAL=${#all_dirs[@]}
COUNT=0
FAIL_COUNT=0

log "Found $TOTAL folder(s) to process across ${#root_dirs[@]} root(s)."

START_TIME=$(date +%s)

for dir in "${all_dirs[@]}"; do
    COUNT=$((COUNT+1))
    log "PROGRESS: [$COUNT/$TOTAL] Processing $dir"
    if ! run_block "$dir"; then
        FAIL_COUNT=$((FAIL_COUNT+1))
    fi
done

END_TIME=$(date +%s)
log "ALL DONE"
log "TOTAL FOLDERS: $TOTAL | FAILED: $FAIL_COUNT"
log "TOTAL TIME: $((END_TIME-START_TIME)) seconds"

if [ "$FAIL_COUNT" -gt 0 ]; then
    exit 1
fi
