#!/bin/bash
#SBATCH --job-name=xfer_4d_1m_eORCA025
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_4d_1m_to_zarr.slurm
#
# Description: SLURM script to send the eORCA025 monthly-mean 
# 4-D output variables to the zarr.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2026-03-21
#
# ==============================================================
# -- Input arguments to OceanDataStore -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 monthly mean output files:
filedir=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA025_JRA55/exp_npd_v1
filepath_gridT=$filedir/eORCA025_1m_grid_T_*.nc
filepath_gridU=$filedir/eORCA025_1m_grid_U_*.nc
filepath_gridV=$filedir/eORCA025_1m_grid_V_*.nc
filepath_gridW=$filedir/eORCA025_1m_grid_W_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=./dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-jra55v1

# -- Python Environment -- #
# Activate Python environment:
source /dssgfs01/working/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA025 monthly mean outputs to JASMIN OS -- #
variables_gridT="time_centered_bounds time_counter_bounds e3t thetao_con so_abs rsdo"
echo "In Progress: Sending eORCA025 4-D T1m variables to zarr..."
ods send_to_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1m_4d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                 -cs '{"x":720, "y":603, "deptht":25}' -v $variables_gridT  -dc $dask_config_json \
                 -zv 3
                    

variables_gridU="time_centered_bounds time_counter_bounds e3u uo uo_eiv u2o umo"
echo "In Progress: Sending eORCA025 4-D U1m variables to zarr..."
ods send_to_zarr -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1m_4d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                 -cs '{"x":720, "y":603, "depthu":25}' -v $variables_gridU  -dc $dask_config_json \
                 -zv 3

variables_gridV="time_centered_bounds time_counter_bounds e3v vo vo_eiv v2o vmo"
echo "In Progress: Sending eORCA025 4-D V1m variables to zarr..."
ods send_to_zarr -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1m_4d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                 -cs '{"x":720, "y":603, "depthv":25}' -v $variables_gridV  -dc $dask_config_json \
                 -zv 3

echo "In Progress: Sending eORCA025 4-D W1m variables to zarr..."
ods send_to_zarr -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1m_4d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                 -cs '{"x":720, "y":603, "depthw":25}'  -dc $dask_config_json \
                 -zv 3

echo "Completed: Sent eORCA025 JRA55-do v1 v1 4-D monthly mean fields to zarr."
