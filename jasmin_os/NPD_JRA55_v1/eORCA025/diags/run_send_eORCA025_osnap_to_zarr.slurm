#!/bin/bash
#SBATCH --job-name=xfer_osnap_eORCA025
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_osnap_to_zarr.slurm
#
# Description: SLURM script to send the eORCA025 section 
# variables for the OSNAP arrays to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2026-03-22
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 OSNAP .zarr store:
filepath_osnap=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA025_JRA55/sections_npd_v1/osnap_npd_v1/eORCA025_1976_2024_OSNAP_section.zarr

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=./dask_config.json

# Define bucket:
bucket=npd-eorca025-jra55v1

# -- Python Environment -- #
# Activate Python environment:
source /dssgfs01/working/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send METRIC eORCA025 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA025 monthly OSNAP outputs to JASMIN object store..."
ods send_to_zarr -f $filepath_osnap -c $store_credentials_json -b $bucket -p M1m/osnap \
                 -dc $dask_config_json -zv 3

echo "Completed: Sent eORCA025 monthly OSNAP outputs to JASMIN object store."
