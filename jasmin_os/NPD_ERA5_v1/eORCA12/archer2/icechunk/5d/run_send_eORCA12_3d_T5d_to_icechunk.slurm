#!/bin/bash
#SBATCH --job-name=xfer_3d_T5d_eORCA12
#SBATCH --time=05:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_3d_T5d_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# 3-D T-point variables to Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-07-18
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

# Define 3-dimensional variables:
variables_gridT="time_centered_bounds time_counter_bounds zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs"

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
# Define the object prefix for the Icechunk:
object_prefix=T5d_3d

# Iterating over years:
for yr in {1980..2024}
do
    # Define the filepath to the 5-day mean T-point files:
    filepath_gridT=$filedir/${yr}/eORCA12_5d_grid_T_${yr}*.nc

    echo "In Progress: Sending eORCA12 ERA-5 v1 T5d 3D variables to Icechunk..."
    ods update_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                        -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                        -cs '{"x":4320, "y":3605}' -v $variables_gridT -dc $dask_config_json \
                        -br "main" -cm "Added eORCA12 ERA5v1 T5d 1976-2024" || { echo "Error: ods update_icechunk failed."; exit 1; }

    echo "Completed: Sent 1976-2024 eORCA12 ERA-5 v1 T5d 3D variables to Icechunk."

done
