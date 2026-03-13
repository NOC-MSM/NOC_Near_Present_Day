#!/bin/bash
#SBATCH --job-name=xfer_1d_S5d_eORCA12
#SBATCH --time=01:00:00
#SBATCH --ntasks=6
#SBATCH --mem=12GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_1d_S5d_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# 1-D scalar variables to Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-07-18
#
# ==============================================================
set -euo pipefail

# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Filepath to JASMIN OS credentials:
store_credentials_json=./credentials/jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=./dask_config.json

# Define the bucket & object prefix for the Icechunk store:
bucket=npd-eorca12-era5v1
object_prefix=S5d_1d

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
# Define the filepath to the 5-day mean scalar files:
filepath_scalar=$filedir/????/eORCA12_5d_scalar_*.nc

echo "In Progress: Sending eORCA12 ERA5v1 S5d variables to Icechunk store..."
ods update_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p $object_prefix \
                    -dc $dask_config_json -br "main" -cm "Added eORCA12 ERA5v1 S5d (1990-01 - 2024-12)"

echo "Completed: Sent 1990-2024 eORCA12 ERA5v1 S5d fields to Icechunk store."
