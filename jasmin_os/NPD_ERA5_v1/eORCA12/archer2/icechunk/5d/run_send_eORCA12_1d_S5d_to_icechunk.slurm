#!/bin/bash
#SBATCH --job-name=xfer_1d_S5d_eORCA12
#SBATCH --time=05:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

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
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Filepath to JASMIN OS credentials:
store_credentials_json=/work/n01/n01/otooth/NOC-NPD/credentials/jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
# Define the filepath to the 5-day mean scalar files:
filepath_scalar=$filedir/????/eORCA12_5d_scalar_*.nc
# Define the object prefix for the Icechunk store:
object_prefix=S5d_1d

echo "In Progress: Sending eORCA12 ERA-5 v1 S5d variables to Icechunk store..."
ods send_to_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p $object_prefix \
                     -br "main" -cm "Added eORCA12 ERA5v1 S5d 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent 1976-2024 eORCA12 ERA-5 v1 S5d fields to Icechunk store."
