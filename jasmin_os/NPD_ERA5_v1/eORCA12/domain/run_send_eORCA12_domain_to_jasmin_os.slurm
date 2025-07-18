#!/bin/bash
#SBATCH --job-name=xfer_domain_cfg_eORCA12
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_domain_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 model domain 
# variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-07-08
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain_cfg=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg_era5v1.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1
# Define the object prefix for the Icechunk store:
object_prefix=domain/domain_cfg

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send ancillary files to JASMIN OS -- #
echo "In Progress: Sending eORCA12 ERA-5 v1 domain_cfg file to Icechunk..."

ods send_to_icechunk -f $filepath_domain_cfg -c $store_credentials_json -b $bucket -p $object_prefix \
                     -cs '{"x":540, "y":721, "nav_lev":25}' \
                     -br "main" -cm "Added eORCA12 ERA5v1 domain_cfg" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 ERA-5 v1 ancillary files to Icechunk."