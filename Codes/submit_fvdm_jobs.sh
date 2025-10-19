#!/bin/bash

mkdir -p logs  # Ensure logs directory exists

# List of density values to run
densities=(2 6 10 14 18 22 26 30 34 38 42 46 50 54 58 62 66 70 74 78 82 86 90 94 98 102 106 110 114 118 122 126 130 134 138 142)
# densities=(42)

for dens in "${densities[@]}"; do
    sbatch unit_fvdm.sh "$dens"
    echo "Submitted job for density $dens"
done