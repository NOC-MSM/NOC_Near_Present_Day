#!/bin/bash
#SBATCH --job-name=xfer_osnap_eORCA025
#SBATCH --partition=test
#SBATCH --time=00:20:00
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
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 OSNAP .nc file:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1/osnap_ERA5_v1
filepath_osnap=$filedir/eORCA025_ERA5_1m_osnap_1976_2024.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# Variables - save entire dataset to single zarr store - and chunks:
variables=consolidated

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send METRIC eORCA025 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA025 monthly OSNAP outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_osnap -c $store_credentials_json -b $bucket -p M1m/osnap -v $variables \
                      -dc $dask_config_json -zv 3 \
                      || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "Completed: Sent eORCA025 monthly OSNAP outputs to JASMIN object store."
