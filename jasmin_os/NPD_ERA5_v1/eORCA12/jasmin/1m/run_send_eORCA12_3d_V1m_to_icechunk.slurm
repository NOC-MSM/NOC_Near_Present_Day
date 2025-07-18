#!/bin/bash
#SBATCH --job-name=xfer_3d_V1m_eORCA12
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_3d_V1m_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# 3-D V-point variables to the Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 monthly mean output files:
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

# Define 3-dimensional variables:
variables_gridV="time_centered_bounds time_counter_bounds tauvo hfy hfy_adv hfy_diff somesatr"

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean V-point files:
filepath_gridV=$filedir/????/eORCA12_1m_grid_V_197*.nc
# Define the object prefix for the Icechunk:
object_prefix=V1m_3d

echo "In Progress: Sending eORCA12 ERA-5 v1 V1m 3D variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p $object_prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                     -cs '{"x":4320, "y":3605}' -v $variables_gridV \
                     -br "main" -cm "Added eORCA12 ERA5v1 V1m 1976-1996" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 ERA-5 v1 V1m 3D variables to Icechunk."
