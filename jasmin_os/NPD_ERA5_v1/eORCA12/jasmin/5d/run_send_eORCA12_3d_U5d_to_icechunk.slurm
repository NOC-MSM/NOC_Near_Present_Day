#!/bin/bash
#SBATCH --job-name=xfer_3d_U5d_eORCA12
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_3d_U5d_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# 3-D U-point variables to Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-07-18
#
# ==============================================================
set -euo pipefail

# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 5-day mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Filepath to JASMIN OS credentials:
store_credentials_json=./credentials/jasmin_object_store_credentials.json

# Define bucket and object prefix:
bucket=npd-eorca12-era5v1
object_prefix=U5d_3d

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# Define 3-dimensional variables:
variables_gridU="time_centered_bounds time_counter_bounds tauuo umo_vint hfx hfx_adv hfx_diff sozosatr"

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
# Define the filepath to the 5-day mean U-point files:
filepath_gridU=$filedir/1990/eORCA12_5d_grid_U_1990*.nc

echo "In Progress: Sending eORCA12 U5d 3-D variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p $object_prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                     -cs '{"x":4320, "y":3605}' -v $variables_gridU -dc $dask_config_json -vc \
                     -br "main" -cm "Added eORCA12 ERA5v1 U5d (1990-01 - 2024-12)"

echo "Completed: Sent 1990 eORCA12 ERA5v1 U5d 3-D variables to Icechunk."

# Iterating over years:
for yr in {1991..2024}
do
    # Define the filepath to the 5-day mean U-point files:
    filepath_gridU=$filedir/${yr}/eORCA12_5d_grid_U_${yr}*.nc

    echo "In Progress: Appending eORCA12 ERA5v1 U5d 3-D ${yr} variables to Icechunk..."
    ods update_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p $object_prefix \
                        -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                        -cs '{"x":4320, "y":3605}' -v $variables_gridU -dc $dask_config_json \
                        -br "main" -cm "Added eORCA12 ERA5v1 U5d (1990-01 - 2024-12)"

    echo "Completed: Appended eORCA12 ERA5v1 U5d 3-D ${yr} variables to Icechunk."

done
