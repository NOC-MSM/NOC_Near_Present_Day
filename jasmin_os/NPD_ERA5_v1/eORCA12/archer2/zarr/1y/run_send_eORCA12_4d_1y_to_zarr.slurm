#!/bin/bash
#SBATCH --job-name=xfer_1y_4d_eORCA12
#SBATCH --time=10:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_4d_1y_to_zarr.slurm
#
# Description: SLURM script to send the eORCA12 monthly-mean 
# output variables to the zarr store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc

# Filepath to eORCA12 annual mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Define filepaths:
filepath_gridT=$filedir/????/eORCA12_1y_grid_T_*.nc
filepath_gridU=$filedir/????/eORCA12_1y_grid_U_*.nc
filepath_gridV=$filedir/????/eORCA12_1y_grid_V_*.nc
filepath_gridW=$filedir/????/eORCA12_1y_grid_W_*.nc

# Define 4-dimensional variables:
variables_gridT="e3t thetao_con so_abs rsdo vohflisf vohfcisf vowflisf"
variables_gridU="e3u uo uo_eiv u2o umo"
variables_gridV="e3v vo vo_eiv v2o vmo"
variables_gridW="depthw_bounds time_centered_bounds time_counter_bounds e3w difvho difvso difvmo avt_evd wmo wo w2o obvfsq"

# Filepath to JASMIN OS credentials:
store_credentials_json=/work/n01/n01/otooth/NOC-NPD/credentials/jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 annual mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA12 T1y 4-D variables to zarr store..."
ods send_to_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":540, "y":721, "deptht":25}' -v $variables_gridT \
                      -dc $dask_config_json -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "Completed: Sent eORCA12 T1y 4-D variables to zarr store."

echo "In Progress: Sending eORCA12 U1y 4-D variables to zarr store..."
ods send_to_zarr -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                      -cs '{"x":540, "y":721, "depthu":25}' -v $variables_gridU \
                      -dc $dask_config_json -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "Completed: Sent eORCA12 U1y 4-D variables to zarr store."

echo "In Progress: Sending eORCA12 V1y 4-D variables to zarr store..."
ods send_to_zarr -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                      -cs '{"x":540, "y":721, "depthv":25}' -v $variables_gridV \
                      -dc $dask_config_json -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "Completed: Sent eORCA12 V1y 4-D variables to zarr store."

echo "In Progress: Sending eORCA12 W1y 4-D variables to zarr store..."
ods send_to_zarr -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":540, "y":721, "depthw":25}' -v $variables_gridW \
                      -dc $dask_config_json -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "Completed: Sent eORCA12 W1y 4-D variables to zarr store."

echo "Completed: Sent eORCA12 annual mean fields to zarr store."
