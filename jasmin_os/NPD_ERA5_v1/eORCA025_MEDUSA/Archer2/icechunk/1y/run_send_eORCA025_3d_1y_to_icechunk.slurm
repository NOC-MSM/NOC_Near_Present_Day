#!/bin/bash
#SBATCH --job-name=xfer_3d_1y_eORCA025
#SBATCH --time=06:00:00
#SBATCH --ntasks=22
#SBATCH --mem=65GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA025_3d_1y_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA025-MEDUSA ERA5v1
# annual-mean 3-D output variables to the Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2026-03-23
#
# ==============================================================
set -euo pipefail

# -- Input arguments to OceanDataStore -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/work/n01/n01/julpal/NEMO/NPD_cray_newer/nemo/cfgs/GLOBAL_QCO_MEDUSA/eORCA025/INPUT/domcfg.nc

# Filepath to eORCA025 annual mean output files:
filedir=/work/n01/n01/julpal/NEMO/NPD_cray_newer/nemo/cfgs/GLOBAL_QCO_MEDUSA/eORCA025/OUTPUT
filepath_gridT=$filedir/eORCA025_1y_grid_T_*.nc
filepath_gridU=$filedir/eORCA025_1y_grid_U_*.nc
filepath_gridV=$filedir/eORCA025_1y_grid_V_*.nc
filepath_gridW=$filedir/eORCA025_1y_grid_W_*.nc
filepath_icemod=$filedir/eORCA025_1y_icemod_*.nc
filepath_scalar=$filedir/eORCA025_1y_scalar_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=./credentials/jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=./dask_config.json

# Define bucket:
bucket=npd-eorca025-medusa-era5v1

# -- Python Environment -- #
# Activate Python environment:
source /home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# === Send eORCA025-MEDUSA ERA5v1 3-D T1y variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA 3-D T1y variables to Icechunk..."

# Define 3-dimensional T1y variables:
variables_gridT="time_centered_bounds time_counter_bounds zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs ttrd_qns_li"

ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y_3d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":1440, "y":1206}' -v $variables_gridT  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 T1y (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 3-D T1y monthly mean fields to Icechunk."
                    
# === Send eORCA025-MEDUSA ERA5v1 3-D U1y variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA 3-D U1y variables to Icechunk..."

# Define 3-dimensional U1y variables:
variables_gridU="time_centered_bounds time_counter_bounds tauuo umo_vint hfx hfx_adv hfx_diff sozosatr"

ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y_3d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                     -cs '{"x":1440, "y":1206}' -v $variables_gridU  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 U1y (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 3-D U1y monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 3-D V1y variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA 3-D V1y variables to Icechunk..."

# Define 3-dimensional V1y variables:
variables_gridV="time_centered_bounds time_counter_bounds tauvo hfy hfy_adv hfy_diff somesatr"

ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y_3d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                     -cs '{"x":1440, "y":1206}' -v $variables_gridV  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 V1y (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 3-D V1y monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 3-D I1y variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA 3-D I1y variables to Icechunk..."

ods send_to_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y_3d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":1440, "y":1206}'  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 I1y (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 3-D I1y monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 3-D S1y variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA 1-D S1y variables to Icechunk..."

ods send_to_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y_1d \
                     -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 S1y (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 1-D S1y monthly mean fields to Icechunk."

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 annual mean fields to Icechunk."
