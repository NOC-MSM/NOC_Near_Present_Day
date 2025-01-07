#!/bin/bash
#SBATCH --job-name=xfer_metric_eORCA025
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_metric_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA025 METRIC section 
# variables for the RAPID, MOVE and SAMBA arrays to the JASMIN
# object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-01-02
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filepath_rapid=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA025_JRA55/metric_npd_v1/eORCA025_1976-01-2024-01_natl_meridional_transports_at_26N.nc
filepath_move=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA025_JRA55/metric_npd_v1/eORCA025_1976-01-2024-01_natl_meridional_transports_at_16N.nc
filepath_samba=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA025_JRA55/metric_npd_v1/eORCA025_1976-01-2024-01_natl_meridional_transports_at_34_5S.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Bucket and object prefix:
bucket=npd-eorca025-jra55v1

# Variables - save entire dataset to single zarr store - and chunks:
variables=compact

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send METRIC eORCA025 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA025 monthly METRIC RAPID 26.5N outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_rapid -c $store_credentials_json -b $bucket -p M1m/rapid_26N -v $variables \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/"}' \
                      -dlc '{"n_workers":10,"threads_per_worker":4,"memory_limit":"3GB"}'

echo "In Progress: Sending eORCA025 monthly METRIC MOVE 16N outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_move -c $store_credentials_json -b $bucket -p M1m/move_16N -v $variables \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/"}' \
                      -dlc '{"n_workers":10,"threads_per_worker":4,"memory_limit":"3GB"}'

echo "In Progress: Sending eORCA025 monthly METRIC SAMBA 34.5S outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_samba -c $store_credentials_json -b $bucket -p M1m/samba_34_5S -v $variables \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/"}' \
                      -dlc '{"n_workers":10,"threads_per_worker":4,"memory_limit":"3GB"}'

echo "Completed: Sent eORCA025 monthly METRIC outputs to JASMIN object store."
