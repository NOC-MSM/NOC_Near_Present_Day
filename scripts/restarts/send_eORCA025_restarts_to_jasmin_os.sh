#!/bin/bash

# ==============================================================
# send_eORCA025_restarts_to_ODS.sh
#
# Description: Script to send eORCA025-ERA5v1 NPD restart .nc
# files to the JASMIN Object Store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-03
#
# ==============================================================
# -- Inputs -- #
# Define restart directory and collect unique restart files:
restart_dir=/dssgfs01/scratch/npd/restarts/eORCA025_ERA5_v1/jasmin_os/
# Define destination path in JASMIN object store:
dest=jasmin-os/npd-eorca025-era5v1/restarts

# Extract restart file list:
cd $restart_dir
restart_files=$(find eORCA025_*_restart.nc)

for file in $restart_files; do
    # Send restart file to JASMIN object store:
    echo "In Progress: Sending $file to ${dest}/${file} JASMIN object store..."
    mc put ${file} ${dest}/${file}
done

echo "Completed: Sent all eORCA025 ERA5v1 restart files to JASMIN object store."
