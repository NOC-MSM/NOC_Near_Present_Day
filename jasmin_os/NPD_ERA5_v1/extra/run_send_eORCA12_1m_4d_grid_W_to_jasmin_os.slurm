#!/bin/bash
#SBATCH --job-name=xfer_W1m_eORCA12
#SBATCH --time=18:00:00
#SBATCH --mem=10G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_1m_grid_W_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# W-point variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-12
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Define 4-dimensional variables:
variables_gridW="depthw_bounds time_centered_bounds time_counter_bounds e3w difvho difvso difvmo avt_evd"

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# Define the object prefix for the JASMIN object store:
object_prefix=W1m

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean W-point files:
# filepath_gridW=$filedir/1976/eORCA12_1m_grid_W_1976*.nc

# echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 W1m variables to JASMIN object store..."
# msm_os send -f $filepath_gridW -c $store_credentials_json -b $bucket -p $object_prefix \
#             -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
#             -cs '{"x":540, "y":721, "depthw":25}' \
#             -zv 3 || { echo "Error: msm_os send failed."; exit 1; }

# echo "Completed: Sent 1976 eORCA12 ERA-5 v1 W1m fields to JASMIN object store."

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1988..1996}
do
    # Define the filepath to the monthly mean W-point files:
    filepath_gridW=$filedir/${yr}/eORCA12_1m_grid_W_*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 W1m variables to JASMIN object store..."
    msm_os update -f $filepath_gridW -c $store_credentials_json -b $bucket -p $object_prefix \
                  -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                  -cs '{"x":540, "y":721, "depthw":25}' -v $variables_gridW \
                  -zv 3 || { echo "Error: msm_os update failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 W1m fields to JASMIN object store."
done
