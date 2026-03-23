#!/bin/bash
#SBATCH --job-name=xfer_1d_S1m_eORCA12
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_1d_S1m_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# 1-D scalar variables to the Icechunk store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
set -euo pipefail

# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Filepath to JASMIN OS credentials:
store_credentials_json=./credentials/jasmin_object_store_credentials.json

# Defin bucket and object prefix:
bucket=npd-eorca12-era5v1
object_prefix=S1m_1d

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean scalar files:
filepath_scalar=$filedir/????/eORCA12_1m_scalar_*.nc

echo "In Progress: Sending eORCA12 ERA5v1 S1m variables to Icechunk store..."
ods send_to_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p $object_prefix \
                     -br "main" -cm "Added eORCA12 ERA5v1 S1m 1976-2024" -vc

echo "Completed: Sent eORCA12 ERA5v1 S1m fields to Icechunk store."
