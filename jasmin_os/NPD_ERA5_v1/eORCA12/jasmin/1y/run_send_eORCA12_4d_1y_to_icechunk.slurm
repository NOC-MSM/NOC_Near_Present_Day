#!/bin/bash
#SBATCH --job-name=xfer_4d_1y_eORCA12
#SBATCH --time=20:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_4d_1y_to_zarr.slurm
#
# Description: SLURM script to send the eORCA12 monthly-mean 
# output variables to Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Define filepaths:
filepath_gridT=$filedir/????/eORCA12_1y_grid_T_*.nc
filepath_gridU=$filedir/????/eORCA12_1y_grid_U_*.nc
filepath_gridV=$filedir/????/eORCA12_1y_grid_V_*.nc
filepath_gridW=$filedir/????/eORCA12_1y_grid_W_*.nc

# Define 4-dimensional variables:
variables_gridT="time_centered_bounds time_counter_bounds e3t thetao_con so_abs rsdo vohflisf vohfcisf vowflisf"
variables_gridU="time_centered_bounds time_counter_bounds e3u uo u2o umo"
variables_gridV="time_centered_bounds time_counter_bounds e3v vo v2o vmo"
variables_gridW="depthw_bounds time_centered_bounds time_counter_bounds e3w difvho difvso difvmo avt_evd wmo wo w2o obvfsq"

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 annual mean outputs to JASMIN OS -- #
# echo "In Progress: Sending eORCA12 T1y 4-D variables to Icechunk..."
# ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y_4d \
#                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
#                      -cs '{"x":540, "y":721, "deptht":25}' -v $variables_gridT \
#                      -br "main" -cm "Added eORCA12 ERA5v1 T1y 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

# echo "Completed: Sent eORCA12 T1y 4-D variables to Icechunk."

# echo "In Progress: Sending eORCA12 U1y 4-D variables to Icechunk..."
# ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y_4d \
#                      -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
#                      -cs '{"x":540, "y":721, "depthu":25}' -v $variables_gridU \
#                      -br "main" -cm "Added eORCA12 ERA5v1 U1y 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

# echo "Completed: Sent eORCA12 U1y 4-D variables to Icechunk."

# echo "In Progress: Sending eORCA12 V1y 4-D variables to Icechunk..."
# ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y_4d \
#                      -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
#                      -cs '{"x":540, "y":721, "depthv":25}' -v $variables_gridV \
#                      -br "main" -cm "Added eORCA12 ERA5v1 V1y 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

# echo "Completed: Sent eORCA12 V1y 4-D variables to Icechunk."

echo "In Progress: Sending eORCA12 W1y 4-D variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1y_4d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":540, "y":721, "depthw":25}' -v $variables_gridW \
                     -br "main" -cm "Added eORCA12 ERA5v1 W1y 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 W1y 4-D variables to Icechunk."

echo "Completed: Sent eORCA12 annual mean fields to Icechunk."
