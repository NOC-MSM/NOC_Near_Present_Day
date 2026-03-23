#!/bin/bash
#SBATCH --job-name=xfer_metric_eORCA025
#SBATCH --partition=test
#SBATCH --time=00:20:00
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
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1/metric_ERA5_v1
filepath_rapid=$filedir/eORCA025_ERA5_1976-01-2024-09_natl_meridional_transports_at_26N.nc
filepath_move=$filedir/eORCA025_ERA5_1976-01-2024-09_natl_meridional_transports_at_16N.nc
filepath_samba=$filedir/eORCA025_ERA5_1976-01-2024-09_natl_meridional_transports_at_34_5S.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=..../jasmin_os/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# Variables - save entire dataset to single zarr store - and chunks:
variables=consolidated

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send METRIC eORCA025 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA025 ERA-5 v1 monthly METRIC RAPID 26.5N outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_rapid -c $store_credentials_json -b $bucket -p M1m/rapid_26N \
                      -v $variables -dc $dask_config_json -zv 3 \
                      || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA025 ERA-5 v1 monthly METRIC MOVE 16N outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_move -c $store_credentials_json -b $bucket -p M1m/move_16N \
                      -v $variables -dc $dask_config_json -zv 3 \
                      || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA025 ERA-5 v1 monthly METRIC SAMBA 34.5S outputs to JASMIN object store..."
msm_os send_with_dask -f $filepath_samba -c $store_credentials_json -b $bucket -p M1m/samba_34_5S \
                      -v $variables -dc $dask_config_json -zv 3 \
                      || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "Completed: Sent eORCA025 ERA-5 v1 monthly METRIC outputs to JASMIN object store."
