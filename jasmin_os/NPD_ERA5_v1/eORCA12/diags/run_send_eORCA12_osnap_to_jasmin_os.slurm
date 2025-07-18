#!/bin/bash
#SBATCH --job-name=xfer_osnap_eORCA025
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_osnap_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA025 section 
# variables for the OSNAP arrays to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA025 OSNAP .nc file:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1/osnap_ERA5_v1
filepath_osnap=eORCA025_ERA5_1976_2024_OSNAP_section.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# Variables - save entire dataset to single zarr store - and chunks:
variables=compact

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send METRIC eORCA025 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA025 monthly OSNAP outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_osnap -c $store_credentials_json -b $bucket -p M1m/osnap -v $variables \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/osnap/M1m/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/osnap/M1m/"}' \
                      -dlc '{"n_workers":15,"threads_per_worker":4,"memory_limit":"3GB"}'

echo "Completed: Sent eORCA025 monthly OSNAP outputs to JASMIN object store."
