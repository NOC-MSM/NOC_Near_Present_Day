#!/bin/bash

# ==============================================================
# rebuild_eORCA025_ocean_restarts.slurm
#
# Description: Script to rebuild NPD-eORCA025-ERA5v1 ocean restart
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
restart_dir=$nemo_dir/cfgs/GLOBAL_QCO/eORCA025/RESTARTS/
# Define output directory for rebuilt restart files:
output_dir=/dssgfs01/scratch/npd/restarts/eORCA025_ERA5_v1

# Define number of domains:
ndom=1326
# Define start year:
yr=1976
# Define initial counter value:
n=0

# Get restart file list:
cd $restart_dir
restart_files=$(find eORCA025_ERA5_*_restart_0000.nc)

# Iterate over restart files and rebuild:
for restart_file in $restart_files; do
    # Get prefix for ocean restart file:
    restart_prefix=${restart_file::30}
    
    # Rebuild ocean restart file:
    echo "In Progress: Rebuilding -> $restart_prefix"
    $nemo_dir/tools/REBUILD_NEMO/rebuild_nemo -m $restart_prefix $ndom
    
    # Move rebuilt ocea restart file to NPD restarts directory:
    if [[ $n -lt 2 ]]; then
        # Spin-up restart files:
        mv $restart_prefix.nc $output_dir/${restart_prefix::8}_${yr}12_${restart_prefix:14}.nc
    else
        mv $restart_prefix.nc $output_dir/${restart_prefix::8}_${yr}12_${restart_prefix:14}.nc
        # Increment year:
        yr=$((yr+1))
    fi

    # Increment counter:
    n=$((n+1))
done
