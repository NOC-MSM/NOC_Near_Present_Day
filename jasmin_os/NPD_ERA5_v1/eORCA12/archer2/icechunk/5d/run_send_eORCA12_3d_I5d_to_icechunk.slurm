#!/bin/bash
#SBATCH --job-name=xfer_3d_I5d_eORCA12
#SBATCH --time=08:00:00
#SBATCH --ntasks=12
#SBATCH --mem=36GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_3d_I5d_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# icemod variables to Icechunk.
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
object_prefix=I5d_3d

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
Define the filepath to the 5-day mean icemod files:
filepath_icemod=$filedir/1990/eORCA12_5d_icemod_1990*.nc

echo "In Progress: Sending eORCA12 I5d 3-D variables to Icechunk..."
ods send_to_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":4320, "y":3605}' -dc $dask_config_json -vc \
                     -br "main" -cm "Added eORCA12 ERA5v1 I5d (1990-01 - 2024-12)"

echo "Completed: Sent 1990-2024 eORCA12 ERA5v1 I5d 3D variables to Icechunk."

# -- Update eORCA12 5-day mean outputs in JASMIN OS -- #
for yr in {1991..2024}
do
    # Define the filepath to the 5-day mean icemod files:
    filepath_icemod=$filedir/${yr}/eORCA12_5d_icemod_${yr}*.nc

    echo "In Progress: Appending eORCA12 ERA5v1 I5d 3-D ${yr} variables to Icechunk store..."
    ods update_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                        -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                        -cs '{"x":4320, "y":3605}' -dc $dask_config_json \
                        -br "main" -cm "Added eORCA12 ERA5v1 I5d (1990-01 - 2024-12)"

    echo "Completed: Appended eORCA12 ERA5v1 I5d 3-D ${yr} fields to Icechunk store."

done
