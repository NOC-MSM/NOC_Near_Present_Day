#!/bin/bash
#SBATCH --job-name=xfer_T1m_3d_eORCA12
#SBATCH --time=10:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_3d_T1m_to_zarr.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# 3-D T-point variables to the zarr store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Define 3-dimensional variables:
variables_gridT="zos zossq friver hfto hfsr hfns hfds rlntds hfss hfls hfempds rsntds mlotst mlotstsq mlotstmax mlotstmin ficeberg berg_latent_heat_flux pbo tos_con tossq_con sos_abs sossq_abs evs prsn hfrainds hfevapds hflx_rnf sfdsi fsitherm empmr tnpeo snowpre snow_ai_cea sowaflup sosafldo somixhgt somxl010 somxzint1 soicecov sowindsp sohflisf sohfcisf sowflisf ocontempmint somint_abs sbt_con sbs_abs"

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

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean T-point files:
filepath_gridT=$filedir/1976/eORCA12_1m_grid_T_1976*.nc
# Define the object prefix for the zarr store:
object_prefix=T1m

echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 T1m 3D variables to zarr store..."
ods send_to_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":2160, "y":1803}' -v $variables_gridT \
                      -dc $dask_config_json -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "Completed: Sent 1976 eORCA12 ERA-5 v1 T1m 3D variables to zarr store."

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1977..2002}
do
    # Define the filepath to the monthly mean T-point files:
    filepath_gridT=$filedir/${yr}/eORCA12_1m_grid_T_${yr}*.nc

    echo "In Progress: Sending $yr eORCA12 ERA-5 v1 T1m 3D variables to zarr store..."
    ods update_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                            -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                            -cs '{"x":2160, "y":1803}' -v $variables_gridT \
                            -dc $dask_config_json -zv 3 || { echo "Error: ods update_zarr failed."; exit 1; }

    echo "Completed: Sent $yr eORCA12 ERA-5 v1 T1m 3D variables to zarr store."
done
