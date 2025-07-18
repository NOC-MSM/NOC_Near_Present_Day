#!/bin/bash
#SBATCH --job-name=xfer_T1m_4d_eORCA12
#SBATCH --time=10:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_4d_T1m_to_zarr.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# T-point 4D variables to the zarr store.
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

# Define 4-dimensional variables:
variables_gridT="e3t thetao_con so_abs rsdo vohflisf vohfcisf vowflisf"

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

echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 T1m 4D variables to zarr store..."
ods send_to_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":540, "y":721, "deptht":25}' -v $variables_gridT \
                      -dc $dask_config_json -zv 3 || { echo "Error: ods send_to_zarr failed."; exit 1; }

echo "Completed: Sent 1976 eORCA12 ERA-5 v1 T1m 4D variables to zarr store."

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1996..2002}
do
    # Define the filepath to the monthly mean T-point files:
    filepath_gridT=$filedir/${yr}/eORCA12_1m_grid_T_${yr}*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 T1m 4D variables to zarr store..."
    ods update_zarr -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                            -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                            -cs '{"x":540, "y":721, "deptht":25}' -v $variables_gridT \
                            -dc $dask_config_json -zv 3 || { echo "Error: ods update_zarr failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 T1m 4D variables to zarr store."
done
