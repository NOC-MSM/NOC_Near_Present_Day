#!/bin/bash
#SBATCH --job-name=xfer_4d_W5d_eORCA12
#SBATCH --time=18:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_4d_W5d_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# W-point variables to Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-07-18
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 5-day mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# Define 4-dimensional variables:
variables_gridW="depthw_bounds time_centered_bounds time_counter_bounds e3w difvho difvso difvmo avt_evd wmo wo w2o obvfsq"

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
# Define the filepath to the 5-day mean W-point files:
filepath_gridW=$filedir/????/eORCA12_5d_grid_W_197*.nc
# Define the object prefix for the Icechunk:
object_prefix=W5d_4d

echo "In Progress: Sending eORCA12 W5d 4-D variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridW -c $store_credentials_json -b $bucket -p $object_prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":540, "y":721, "depthw":25}' -v $variables_gridW -vc \
                     -br "main" -cm "Added eORCA12 ERA5v1 W5d 1976-1979" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent 1976-1979 eORCA12 ERA-5 v1 W5d 4D variables to Icechunk."
