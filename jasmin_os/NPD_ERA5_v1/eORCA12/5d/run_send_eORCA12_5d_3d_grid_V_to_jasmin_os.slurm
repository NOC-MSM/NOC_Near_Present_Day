#!/bin/bash
#SBATCH --job-name=xfer_V5d_3d_eORCA12
#SBATCH --time=10:00:00
#SBATCH --ntasks=15
#SBATCH --mem=38GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_1m_3d_grid_V_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# 3-D V-point variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc

# Filepath to eORCA12  5d mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Define 3-dimensional variables:
variables_gridV="depthv_bounds time_centered_bounds time_counter_bounds tauvo somesatr"

# Filepath to JASMIN OS credentials:
store_credentials_json=/work/n01/n01/otooth/NOC-NPD/credentials/jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA12/5d/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Define the filepath to the  5d mean V-point files:
filepath_gridV=$filedir/1976/eORCA12_5d_grid_V_1976*.nc
# Define the object prefix for the JASMIN object store:
object_prefix=V5d

echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 V5d 3D variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridV -c $store_credentials_json -b $bucket -p $object_prefix \
                      -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                      -cs '{"x":2160, "y":1803}' -v $variables_gridV \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "Completed: Sent 1976 eORCA12 ERA-5 v1 V5d 3D variables to JASMIN object store."

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1977..1996}
do
    # Define the filepath to the 5d mean V-point files:
    filepath_gridV=$filedir/${yr}/eORCA12_5d_grid_V_${yr}*.nc

    echo "In Progress: Sending $yr eORCA12 ERA-5 v1 V5d 3D variables to JASMIN object store..."
    msm_os update_with_dask -f $filepath_gridV -c $store_credentials_json -b $bucket -p $object_prefix \
                            -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                            -cs '{"x":2160, "y":1803}' -v $variables_gridV \
                            -dc $dask_config_json -zv 3 || { echo "Error: msm_os update_with_dask failed."; exit 1; }

    echo "Completed: Sent $yr eORCA12 ERA-5 v1 V5d 3D variables to JASMIN object store."
done
