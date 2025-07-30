#!/bin/bash
#SBATCH --job-name=xfer_T1m_eORCA12
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_1m_3d_grid_T_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# 3-D T-point variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Define 3-dimensional variables:
variables_gridT="zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs"

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean T-point files:
filepath_gridT=$filedir/1976/eORCA12_1m_grid_T_1976*.nc
# Define the object prefix for the JASMIN object store:
object_prefix=T1m

echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 T1m 3D variables to JASMIN object store..."
msm_os send -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
            -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
            -cs '{"x":2160, "y":1803}' -v $variables_gridT \
            -zv 3 || { echo "Error: msm_os send failed."; exit 1; }

echo "Completed: Sent 1976 eORCA12 ERA-5 v1 T1m 3D variables to JASMIN object store."

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1991..1996}
do
    # Define the filepath to the monthly mean T-point files:
    filepath_gridT=$filedir/${yr}/eORCA12_1m_grid_T_${yr}*.nc

    echo "In Progress: Sending $yr eORCA12 ERA-5 v1 T1m 3D variables to JASMIN object store..."
    msm_os update -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                  -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                  -cs '{"x":2160, "y":1803}' -v $variables_gridT \
                  -zv 3 || { echo "Error: msm_os update failed."; exit 1; }

    echo "Completed: Sent $yr eORCA12 ERA-5 v1 T1m 3D variables to JASMIN object store."
done
