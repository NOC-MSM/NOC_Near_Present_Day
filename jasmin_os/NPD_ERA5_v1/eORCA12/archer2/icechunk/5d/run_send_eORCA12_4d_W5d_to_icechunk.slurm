#!/bin/bash
#SBATCH --job-name=xfer_4d_W5d_eORCA12
#SBATCH --time=10:00:00
#SBATCH --ntasks=25
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

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

# Define bucket and object prefix:
bucket=npd-eorca12-era5v1
object_prefix=W5d_4d

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# Define 4-dimensional variables:
variables_gridW="time_centered_bounds time_counter_bounds depthw_bounds e3w wo wmo obvfsq"

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
# Define the filepath to the 5-day mean W-point files:
filepath_gridW=$filedir/1990/eORCA12_5d_grid_W_1990*.nc

echo "In Progress: Sending eORCA12 W5d 4-D variables to Icechunk..."
ods update_icechunk -f $filepath_gridW -c $store_credentials_json -b $bucket -p $object_prefix \
                    -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                    -cs '{"x":540, "y":721, "depthw":25}' -v $variables_gridW -dc $dask_config_json -vc \
                    -br "main" -cm "Added eORCA12 ERA5v1 W5d (1990-01 - 2024-12)"

echo "Completed: Sent 1990 eORCA12 ERA-5 v1 W5d 4-D variables to Icechunk."

# Iterating over years:
for yr in {1991..2024}
do
    # Define the filepath to the 5-day mean W-point files:
    filepath_gridW=$filedir/${yr}/eORCA12_5d_grid_W_${yr}*.nc

    echo "In Progress: Appending eORCA12 W5d 4-D ${yr} variables to Icechunk..."
    ods update_icechunk -f $filepath_gridW -c $store_credentials_json -b $bucket -p $object_prefix \
                        -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                        -cs '{"x":540, "y":721, "depthw":25}' -v $variables_gridW -dc $dask_config_json \
                        -br "main" -cm "Added eORCA12 ERA5v1 W5d (1990-01 - 2024-12)"

    echo "Completed: Appended eORCA12 ERA-5 v1 W5d 4-D ${yr} variables to Icechunk."

done
