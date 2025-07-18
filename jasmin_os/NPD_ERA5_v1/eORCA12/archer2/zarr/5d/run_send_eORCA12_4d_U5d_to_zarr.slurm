#!/bin/bash
#SBATCH --job-name=xfer_U5d_4d_eORCA12
#SBATCH --time=12:00:00
#SBATCH --ntasks=10
#SBATCH --mem=25GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_4d_U5d_to_zarr.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# U-point 4D variables to the zarr store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc

# Filepath to eORCA12 5d mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Define 4-dimensional variables:
variables_gridU="e3u uo uo_eiv u2o umo"

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

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Define the filepath to the 5d mean U-point files:
filepath_gridU=$filedir/1977/eORCA12_5d_grid_U_1977*.nc
# Define the object prefix for the zarr store:
object_prefix=U5d

echo "In Progress: Sending 1977 eORCA12 ERA-5 v1 U5d 4D variables to zarr store..."
ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p $object_prefix \
                      -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                      -cs '{"x":540, "y":721, "depthu":25}' -v $variables_gridU \
                      -dc $dask_config_json -zv 3 || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent 1977 eORCA12 ERA-5 v1 U5d 4D variables to zarr store."

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1978..1996}
do
    # Define the filepath to the 5d mean U-point files:
    filepath_gridU=$filedir/${yr}/eORCA12_5d_grid_U_${yr}*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 U5d 4D variables to zarr store..."
    ods update_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p $object_prefix \
                            -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                            -cs '{"x":540, "y":721, "depthu":25}' -v $variables_gridU \
                            -dc $dask_config_json -zv 3 || { echo "Error: ods update_icechunk failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 U5d 4D variables to zarr store."
done
