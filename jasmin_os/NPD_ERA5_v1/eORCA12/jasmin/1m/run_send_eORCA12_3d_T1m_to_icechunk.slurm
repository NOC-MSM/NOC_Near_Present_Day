#!/bin/bash
#SBATCH --job-name=xfer_3d_T1m_eORCA12
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_3d_T1m_to_Icechunk.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# 3-D T-point variables to the Icechunk.
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

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# Define 3-dimensional variables:
variables_gridT="time_centered_bounds time_counter_bounds zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs"

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean T-point files:
filepath_gridT=$filedir/????/eORCA12_1m_grid_T_197*.nc
# Define the object prefix for the Icechunk:
object_prefix=T1m_3d

echo "In Progress: Sending eORCA12 ERA-5 v1 T1m 3D variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":4320, "y":3605}' -v $variables_gridT \
                     -br "main" -cm "Added eORCA12 ERA5v1 T1m 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA12 ERA-5 v1 T1m 3D variables to Icechunk."
