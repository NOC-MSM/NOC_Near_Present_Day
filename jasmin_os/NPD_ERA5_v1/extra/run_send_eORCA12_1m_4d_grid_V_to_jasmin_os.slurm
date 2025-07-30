#!/bin/bash
#SBATCH --job-name=xfer_V1m_eORCA12
#SBATCH --time=18:00:00
#SBATCH --mem=10G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_1m_4d_grid_V_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# V-point 4D variables to the JASMIN object store.
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
variables_gridV="e3v vo vo_eiv v2o vmo"

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
# Define the filepath to the monthly mean V-point files:
# filepath_gridV=$filedir/1976/eORCA12_1m_grid_V_1976*.nc
# Define the object prefix for the JASMIN object store:
object_prefix=V1m

# echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 V1m 4D variables to JASMIN object store..."
# msm_os send -f $filepath_gridV -c $store_credentials_json -b $bucket -p $object_prefix \
#             -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
#             -cs '{"x":540, "y":721, "depthv":25}' -v $variables_gridV \
#             -zv 3 || { echo "Error: msm_os send failed."; exit 1; }

# echo "Completed: Sent 1976 eORCA12 ERA-5 v1 V1m 4D variables to JASMIN object store."

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1995..1996}
do
    # Define the filepath to the monthly mean V-point files:
    filepath_gridV=$filedir/${yr}/eORCA12_1m_grid_V_${yr}*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 V1m 4D variables to JASMIN object store..."
    msm_os update -f $filepath_gridV -c $store_credentials_json -b $bucket -p $object_prefix \
                  -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                  -cs '{"x":540, "y":721, "depthv":25}' -v $variables_gridV \
                  -zv 3 || { echo "Error: msm_os update failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 V1m 4D variables to JASMIN object store."
done
