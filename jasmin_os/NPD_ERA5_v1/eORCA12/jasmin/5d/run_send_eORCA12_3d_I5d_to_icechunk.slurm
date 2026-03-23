#!/bin/bash
#SBATCH --job-name=xfer_3d_I5d_eORCA12
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

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
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 5-day mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Filepath to JASMIN OS credentials:
store_credentials_json=./credentials/jasmin_object_store_credentials.json

# Defin bucket and object prefix:
bucket=npd-eorca12-era5v1
object_prefix=I5d_3d

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
# Define the filepath to the 5-day mean icemod files:
filepath_icemod=$filedir/1990/eORCA12_5d_icemod_1990*.nc

echo "In Progress: Sending eORCA12 ERA5v1 I5d variables to Icechunk store..."
ods send_to_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":4320, "y":3605}' -vc \
                     -br "main" -cm "Added eORCA12 ERA5v1 I5d (1990-01 - 2024-12)"

echo "Completed: Sent eORCA12 ERA5v1 I5d fields to Icechunk store."

echo "In Progress: Appending eORCA12 ERA5v1 I5d variables to Icechunk store..."
for yr in {1991..2024}
do
  # Define the filepath to the 5-day mean icemod files:
  filepath_icemod=$filedir/${yr}/eORCA12_5d_icemod_${yr}*.nc

  echo "In Progress: Appending eORCA12 ERA5v1 I5d ${yr} variables to Icechunk..."
  ods update_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":4320, "y":3605}' \
                      -br "main" -cm "Added eORCA12 ERA5v1 I5d (1990-01 - 2024-12)"

    echo "Completed: Appended eORCA12 ERA5v1 I5d ${yr} fields to Icechunk store."
done
