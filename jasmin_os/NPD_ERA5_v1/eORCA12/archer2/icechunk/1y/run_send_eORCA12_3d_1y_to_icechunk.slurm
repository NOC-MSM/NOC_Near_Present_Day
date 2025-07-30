#!/bin/bash
#SBATCH --job-name=xfer_1y_3d_eORCA12
#SBATCH --time=12:00:00
#SBATCH --ntasks=12
#SBATCH --mem=25GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_3d_1y_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA12 monthly-mean 
# 3-dimensional output variables to Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-18
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA12/domain_cfg_xy_adjusted.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

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

for yr in {1980..2023}
do

# Define filepaths:
filepath_gridT=$filedir/????/eORCA12_1y_grid_T_${yr}*.nc
filepath_gridU=$filedir/????/eORCA12_1y_grid_U_${yr}*.nc
filepath_gridV=$filedir/????/eORCA12_1y_grid_V_${yr}*.nc
filepath_icemod=$filedir/????/eORCA12_1y_icemod_${yr}*.nc
filepath_scalar=$filedir/????/eORCA12_1y_scalar_${yr}*.nc

# Define 3-dimensional variables:
variables_gridT="time_centered_bounds time_counter_bounds zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs"
variables_gridU="time_centered_bounds time_counter_bounds tauuo umo_vint hfx hfx_adv hfx_diff sozosatr"
variables_gridV="time_centered_bounds time_counter_bounds tauvo hfy hfy_adv hfy_diff somesatr"

# -- Send eORCA12 annual mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA12 S1y 1-D variables to Icechunk..."
ods update_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y_1d \
                    -br "main" -cm "Added eORCA12 ERA5v1 S1y 1976-2024" || { echo "Error: ods update_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 S1y 1-D variables to Icechunk."

echo "In Progress: Sending eORCA12 T1y 3-D variables to Icechunk..."
ods update_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y_3d \
                    -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                    -cs '{"x":4320, "y":3605}' -v $variables_gridT \
                    -br "main" -cm "Added eORCA12 ERA5v1 T1y 1976-2024" || { echo "Error: ods update_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 T1y 3-D variables to Icechunk."

echo "In Progress: Sending eORCA12 U1y 3-D variables to Icechunk..."
ods update_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y_3d \
                    -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                    -cs '{"x":4320, "y":3605}' -v $variables_gridU \
                    -br "main" -cm "Added eORCA12 ERA5v1 U1y 1976-2024" || { echo "Error: ods update_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 U1y 3-D variables to Icechunk."

echo "In Progress: Sending eORCA12 V1y 3-D variables to Icechunk..."
ods update_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y_3d \
                    -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                    -cs '{"x":4320, "y":3605}' -v $variables_gridV \
                    -br "main" -cm "Added eORCA12 ERA5v1 V1y 1976-2024" || { echo "Error: ods update_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 V1y 3-D variables to Icechunk."

echo "In Progress: Sending eORCA12 I1y 3-D variables to Icechunk..."
ods update_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y_3d \
                    -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                    -cs '{"x":4320, "y":3605}' \
                    -br "main" -cm "Added eORCA12 ERA5v1 I1y 1976-2024" || { echo "Error: ods update_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 I1y 3-D variables to Icechunk."

done

echo "Completed: Sent eORCA12 annual mean fields to Icechunk."
