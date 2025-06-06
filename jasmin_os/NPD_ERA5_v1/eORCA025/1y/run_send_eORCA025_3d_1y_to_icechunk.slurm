#!/bin/bash
#SBATCH --job-name=xfer_3d_1y_eORCA025
#SBATCH --partition=compute
#SBATCH --time=05:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_3d_1y_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA025 annual-mean 
# 3-D output variables to the Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-05-26
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 annual mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1
filepath_gridT=$filedir/eORCA025_ERA5_1y_grid_T_*.nc
filepath_gridU=$filedir/eORCA025_ERA5_1y_grid_U_*.nc
filepath_gridV=$filedir/eORCA025_ERA5_1y_grid_V_*.nc
filepath_gridW=$filedir/eORCA025_ERA5_1y_grid_W_*.nc
filepath_icemod=$filedir/eORCA025_ERA5_1y_icemod_*.nc
filepath_scalar=$filedir/eORCA025_ERA5_1y_scalar_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=/dssgfs01/working/otooth/AtlantiS/jasmin_os/credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=/dssgfs01/scratch/otooth/NOC_NPD/20250512_eORCA1_ERA5v1/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA025/1y/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 annual mean outputs to JASMIN OS -- #
variables_gridT="time_centered_bounds time_counter_bounds zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs ttrd_qns_li"

echo "In Progress: Sending eORCA025 3-D T1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y_3d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":1440, "y":1206}' -v $variables_gridT  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 T1y 1976-2024" -vc || { echo "Error: ods update_icechunk failed."; exit 1; }
                    

variables_gridU="time_centered_bounds time_counter_bounds tauuo umo_vint hfx hfx_adv hfx_diff sozosatr"

echo "In Progress: Sending eORCA025 3-D U1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y_3d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                     -cs '{"x":1440, "y":1206}' -v $variables_gridU  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 U1y 1976-2024" -vc || { echo "Error: ods update_icechunk failed."; exit 1; }

variables_gridV="time_centered_bounds time_counter_bounds tauvo hfy hfy_adv hfy_diff somesatr"

echo "In Progress: Sending eORCA025 3-D V1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y_3d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                     -cs '{"x":1440, "y":1206}' -v $variables_gridV  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 V1y 1976-2024" -vc || { echo "Error: ods update_icechunk failed."; exit 1; }

# echo "In Progress: Sending eORCA025 I1y variables to Icechunk..."
# ods send_to_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y \
#                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
#                      -cs '{"x":1440, "y":1206}'  -dc $dask_config_json \
#                      -br "main" -cm "Added eORCA025 ERA5v1 I1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

# echo "In Progress: Sending eORCA025 S1y variables to Icechunk..."
# ods send_to_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y \
#                      -dc $dask_config_json \
#                      -br "main" -cm "Added eORCA025 ERA5v1 S1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

# echo "Completed: Sent eORCA025 ERA-5 v1 annual mean fields to Icechunk."
