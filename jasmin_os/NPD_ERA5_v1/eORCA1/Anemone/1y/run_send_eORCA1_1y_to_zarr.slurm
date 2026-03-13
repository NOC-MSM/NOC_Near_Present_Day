#!/bin/bash
#SBATCH --job-name=xfer_1y_eORCA1
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_1y_to_zarr.slurm
#
# Description: SLURM script to send the eORCA1 ERA-5 v1 annual-mean 
# output variables to the Zarr.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA1 ancillary file:
filepath_domain=/dssgfs01/scratch/npd/simulations/Domains/eORCA1/domain_cfg.nc

# Filepath to eORCA1 annual mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA1_ERA5_v1
filepath_gridT=$filedir/eORCA1_ERA5_1y_grid_T_*.nc
filepath_gridU=$filedir/eORCA1_ERA5_1y_grid_U_*.nc
filepath_gridV=$filedir/eORCA1_ERA5_1y_grid_V_*.nc
filepath_gridW=$filedir/eORCA1_ERA5_1y_grid_W_*.nc
filepath_icemod=$filedir/eORCA1_ERA5_1y_icemod_*.nc
filepath_scalar=$filedir/eORCA1_ERA5_1y_scalar_*.nc

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

# -- Send eORCA1 ERA-5 v1 annual mean outputs to JASMIN OS -- # 
echo "In Progress: Sending eORCA1 ERA-5 v1 T1y variables to Zarr..."
ods send_to_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y \
                 -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                 -cs '{"x":360,"y":331,"deptht":25}' \
                 -dc $dask_config_json -vs -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 U1y variables to Zarr..."
ods send_to_zarr -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y \
                 -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                 -cs '{"x":360,"y":331,"depthu":25}' \
                 -dc $dask_config_json -vs -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 V1y variables to Zarr..."
ods send_to_zarr -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y \
                 -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                 -cs '{"x":360,"y":331,"depthv":25}' \
                 -dc $dask_config_json -vs -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 W1y variables to Zarr..."
ods send_to_zarr -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1y \
                 -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                 -cs '{"x":360,"y":331,"depthw":25}' \
                 -dc $dask_config_json -vs -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 I1y variables to Zarr..."
ods send_to_zarr -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y \
                 -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                 -cs '{"x":360,"y":331}' \
                 -dc $dask_config_json -vs -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 S1y variables to Zarr..."
ods send_to_zarr -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y \
                 -dc $dask_config_json -vs -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "Completed: Sent eORCA1 ERA-5 v1 annual mean fields to Zarr stores."
