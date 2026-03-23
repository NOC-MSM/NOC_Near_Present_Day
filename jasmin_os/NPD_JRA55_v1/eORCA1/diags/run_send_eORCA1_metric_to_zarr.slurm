#!/bin/bash
#SBATCH --job-name=xfer_metric_eORCA1
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_metric_to_zarr.slurm
#
# Description: SLURM script to send the eORCA1 METRIC section 
# variables for the RAPID, MOVE and SAMBA arrays to the JASMIN
# object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-23
#
# ==============================================================
# -- Input arguments to OceanDataStore -- #
# Filepath to eORCA1 ancillary file:
filepath_rapid=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA1_JRA55/metric_npd_v1/eORCA1_1976-01-2024-01_natl_meridional_transports_at_26N.nc
filepath_move=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA1_JRA55/metric_npd_v1/eORCA1_1976-01-2024-01_natl_meridional_transports_at_16N.nc
filepath_samba=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA1_JRA55/metric_npd_v1/eORCA1_1976-01-2024-01_natl_meridional_transports_at_34_5S.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=./dask_config.json

# Define bucket:
bucket=npd-eorca1-jra55v1

# -- Python Environment -- #
# Activate Python environment:
source /dssgfs01/working/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send METRIC eORCA1 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA1 monthly METRIC RAPID 26.5N outputs to JASMIN object store..."
ods send_to_zarr -f $filepath_rapid -c $store_credentials_json -b $bucket -p M1m/rapid_26N \
                 -dc $dask_config_json -zv 3

echo "In Progress: Sending eORCA1 monthly METRIC MOVE 16N outputs to JASMIN object store..."
ods send_to_zarr -f $filepath_move -c $store_credentials_json -b $bucket -p M1m/move_16N \
                 -dc $dask_config_json -zv 3

echo "In Progress: Sending eORCA1 monthly METRIC SAMBA 34.5S outputs to JASMIN object store..."
ods send_to_zarr -f $filepath_samba -c $store_credentials_json -b $bucket -p M1m/samba_34_5S \
                 -dc $dask_config_json -zv 3

echo "Completed: Sent eORCA1 monthly METRIC outputs to JASMIN object store."
