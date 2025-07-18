#!/bin/bash
#SBATCH --job-name=xfer_3d_I1m_eORCA12
#SBATCH --time=15:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_3d_I1m_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# icemod variables to the Icechunk store.
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

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean icemod files:
filepath_icemod=$filedir/????/eORCA12_1m_icemod_*.nc
# Define the object prefix for the Icechunk store:
object_prefix=I1m_3d

echo "In Progress: Sending eORCA12 ERA-5 v1 I1m variables to Icechunk store..."
ods send_to_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":4320, "y":3605}' \
                     -br "main" -cm "Added eORCA12 ERA5v1 I1m 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 ERA-5 v1 I1m fields to Icechunk store."
