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
# variables for the RAPID, MOVE and SAMBA arrays to Zarr stores
# in the JASMIN Object Store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepaths to eORCA1 ERA-5 v1 diagnostics:
filedir=/dssgfs01/scratch/npd/simulations/eORCA1_ERA5_v1/metric_ERA5_v1
filepath_rapid=$filedir/eORCA1_ERA5_1976-01-2024-09_natl_meridional_transports_at_26N.nc
filepath_move=$filedir/eORCA1_ERA5_1976-01-2024-09_natl_meridional_transports_at_16N.nc
filepath_samba=$filedir/eORCA1_ERA5_1976-01-2024-09_natl_meridional_transports_at_34_5S.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca1-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send METRIC eORCA1 ERA-5 v1 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA1 ERA-5 v1 monthly METRIC RAPID 26.5N outputs to Zarr..."
msm_os send_with_dask -f $filepath_rapid -c $store_credentials_json -b $bucket -p M1m/RAPID_26N \
                      -dc $dask_config_json -zv 3 \

echo "In Progress: Sending eORCA1 ERA-5 v1 monthly METRIC MOVE 16N outputs to Zarr..."
msm_os send_with_dask -f $filepath_move -c $store_credentials_json -b $bucket -p M1m/MOVE_16N \
                      -dc $dask_config_json -zv \

echo "In Progress: Sending eORCA1 ERA-5 v1 monthly METRIC SAMBA 34.5S outputs to Zarr..."
msm_os send_with_dask -f $filepath_samba -c $store_credentials_json -b $bucket -p M1m/SAMBA_34_5S \
                       -dc $dask_config_json -zv 3 \

echo "Completed: Sent eORCA1 ERA-5 v1 monthly METRIC outputs to Zarr stores."
