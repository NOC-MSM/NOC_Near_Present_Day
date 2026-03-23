#!/bin/bash
#SBATCH --job-name=xfer_3d_1y_eORCA025
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_3d_1y_to_zarr.slurm
#
# Description: SLURM script to send the eORCA025 annual-mean 
# 3-D output variables to the Zarr.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2026-03-21
#
# ==============================================================
# -- Input arguments to OceanDataStore -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 annual mean output files:
filedir=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA025_JRA55/exp_npd_v1
filepath_gridT=$filedir/eORCA025_1y_grid_T_*.nc
filepath_gridU=$filedir/eORCA025_1y_grid_U_*.nc
filepath_gridV=$filedir/eORCA025_1y_grid_V_*.nc
filepath_gridW=$filedir/eORCA025_1y_grid_W_*.nc
filepath_icemod=$filedir/eORCA025_1y_icemod_*.nc
filepath_scalar=$filedir/eORCA025_1y_scalar_*.nc

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

# -- Send eORCA025 annual mean outputs to JASMIN OS -- #
variables_gridT="time_centered_bounds time_counter_bounds zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs ttrd_qns_li"

echo "In Progress: Sending eORCA025 3-D T1y variables to zarr..."
ods send_to_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y_3d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                 -cs '{"x":1440, "y":1206}' -v $variables_gridT  -dc $dask_config_json \
                 -zv 3
                    

variables_gridU="time_centered_bounds time_counter_bounds tauuo umo_vint hfx hfx_adv hfx_diff sozosatr"

echo "In Progress: Sending eORCA025 3-D U1y variables to zarr..."
ods send_to_zarr -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y_3d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                 -cs '{"x":1440, "y":1206}' -v $variables_gridU  -dc $dask_config_json \
                 -zv 3

variables_gridV="time_centered_bounds time_counter_bounds tauvo hfy hfy_adv hfy_diff somesatr"

echo "In Progress: Sending eORCA025 3-D V1y variables to zarr..."
ods send_to_zarr -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y_3d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                 -cs '{"x":1440, "y":1206}' -v $variables_gridV  -dc $dask_config_json \
                 -zv 3

echo "In Progress: Sending eORCA025 I1y variables to zarr..."
ods send_to_zarr -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y_3d \
                 -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                 -cs '{"x":1440, "y":1206}'  -dc $dask_config_json \
                 -zv 3

echo "In Progress: Sending eORCA025 S1y variables to zarr..."
ods send_to_zarr -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y_1d \
                 -dc $dask_config_json -zv 3

echo "Completed: Sent eORCA025 JRA55-do v1 v1 annual mean fields to zarr."
