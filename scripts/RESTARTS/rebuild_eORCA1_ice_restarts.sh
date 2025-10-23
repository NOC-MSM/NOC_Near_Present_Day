#!/bin/bash

# ==============================================================
# rebuild_eORCA1_ice_restarts.slurm
#
# Description: Script to rebuild NPD-eORCA1-ERA5v1 ice restart
# files from subdomains in the RESTART/ directory.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-03
#
# ==============================================================
# -- Load Modules -- #
module purge
module load NEMO/prg-env

# -- Inputs -- #
# Define path to NEMO directory in NOC_Near_Present_Day:
nemo_dir=..../NOC_Near_Present_Day/nemo
# Define restart directory and collect unique restart files:
restart_dir=$nemo_dir/cfgs/GLOBAL_QCO/eORCA1/RESTARTS/
# Define output directory for rebuilt restart files:
output_dir=/dssgfs01/scratch/npd/restarts/eORCA1_ERA5_v1/

# Define number of domains:
ndom=552
# Define start year:
yr=1976
# Define initial counter value:
n=0

# Get restart file list:
cd $restart_dir
restart_files=$(find eORCA1_ERA5_*_restart_ice_0000.nc)

# Iterate over restart files and rebuild:
for restart_file in $restart_files; do
    # Get prefix for restart file:
    restart_prefix=${restart_file::28}
    
    # Rebuild ice restart file:
    echo "In Progress: Rebuilding -> $restart_prefix"
    $nemo_dir/tools/REBUILD_NEMO/rebuild_nemo -m $restart_prefix $ndom
    
    # Move rebuilt ice restart file to NPD restarts directory:
    if [[ $n -lt 2 ]]; then
        # Spin-Up restart files:
        mv $restart_prefix.nc $output_dir/${restart_prefix::11}_${yr}12_${restart_prefix:14}.nc
    else
        mv $restart_prefix.nc $output_dir/${restart_prefix::11}_${yr}12_${restart_prefix:14}.nc
        # Increment year:
        yr=$((yr+1))
    fi

    # Increment counter:
    n=$((n+1))
done
