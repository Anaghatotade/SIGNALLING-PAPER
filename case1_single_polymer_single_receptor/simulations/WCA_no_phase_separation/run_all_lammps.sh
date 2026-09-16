#!/bin/bash
# Run all LAMMPS simulations sequentially for the WCA (no phase separation)
# condition of the single-polymer/single-receptor system.
# Author: Anagha
#
# NOTE: fixed two bugs relative to the original script -- it looped over
# "$folders" (the whole array) instead of "$folder" (the loop variable), and
# used an author-specific absolute base_dir. Both are corrected below so the
# script actually visits every N_* folder and works from any checkout location.

# Base directory where N_* folders are located (relative to this script)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_dir="$SCRIPT_DIR"

# List of folders to process
folders=("N_20" "N_50" "N_100" "N_125" "N_150" "N_200" "N_250" "N_300")

# Path to your LAMMPS executable
LAMMPS_CMD="$HOME/Documents/lammps_new/build/lmp -in lammps_sim.in"

# Loop through each folder and run LAMMPS
for folder in "${folders[@]}"; do

    echo "Running simulation in: $folder"

    cd "$base_dir/$folder" || { echo "Cannot cd into $folder"; continue; }

    if [ -f "lammps_sim.in" ]; then
        echo "Starting LAMMPS in $folder..."
        $LAMMPS_CMD > output.log 2>&1
        echo "Finished simulation in $folder"
    else
        echo " lammps_sim.in not found in $folder — skipping."
    fi

    # Go back to main directory
    cd "$base_dir" || exit
done

echo " All simulations completed!"
