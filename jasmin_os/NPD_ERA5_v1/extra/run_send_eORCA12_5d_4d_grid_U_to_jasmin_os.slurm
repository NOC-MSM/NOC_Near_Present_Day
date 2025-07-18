#!/bin/bash
#SBATCH --job-name=xfer_U5d_eORCA12
#SBATCH --time=08:00:00
#SBATCH --mem=10G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_5d_4d_grid_U_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# U-point 4D variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-11
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Define 4-dimensional variables:
variables_gridU="e3u uo uo_eiv u2o umo"

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Define the filepath to the 5d mean U-point files:
filepath_gridU=$filedir/1976/eORCA12_5d_grid_U_1976*.nc
# Define the object prefix for the JASMIN object store:
object_prefix=U5d

echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 U5d 4D variables to JASMIN object store..."
msm_os send -f $filepath_gridU -c $store_credentials_json -b $bucket -p $object_prefix \
            -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
            -cs '{"x":540, "y":721, "depthu":25}' -v $variables_gridU \
            -zv 3 || { echo "Error: msm_os send failed."; exit 1; }

echo "Completed: Sent 1976 eORCA12 ERA-5 v1 U5d 4D variables to JASMIN object store."

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1977..1996}
do
    # Define the filepath to the 5d mean U-point files:
    filepath_gridU=$filedir/${yr}/eORCA12_5d_grid_U*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 U5d 4D variables to JASMIN object store..."
    msm_os update -f $filepath_gridU -c $store_credentials_json -b $bucket -p $object_prefix \
                  -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                  -cs '{"x":540, "y":721, "depthu":25}' -v $variables_gridU \
                  -zv 3 || { echo "Error: msm_os update failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 U5d 4D variables to JASMIN object store."
done
