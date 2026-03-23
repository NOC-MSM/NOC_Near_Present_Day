#!/bin/bash
#SBATCH --job-name=xfer_domain_eORCA1
#SBATCH --partition=test
#SBATCH --time=00:18:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_domain_to_zarr.slurm
#
# Description: SLURM script to send the eORCA1 model domain 
# variables to Zarr store in the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2026-03-22
#
# ==============================================================
set -euo pipefail

# -- Input arguments to OceanDataStore -- #
# Filepath to eORCA1 ancillary files:
filepath_domain_cfg=/dssgfs01/scratch/npd/simulations/Domains/eORCA1/eORCA1_ERA5v1_domain_cfg_mesh_mask.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=./dask_config.json

# Define Bucket:
bucket=npd-eorca1-jra55v1

# -- Python Environment -- #
# Activate Python environment:
source /dssgfs01/working/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send ancillary files to JASMIN OS -- #
echo "In Progress: Sending eORCA1 JRA55 v1 domain_cfg file to zarr..."
ods send_to_zarr -f $filepath_domain_cfg -c $store_credentials_json -b $bucket -p domain_cfg \
                 -cs '{"x":360,"y":331}' -dc $dask_config_json \
                 -zv 3

echo "Completed: Sent eORCA1 JRA55 ancillary files to zarr."