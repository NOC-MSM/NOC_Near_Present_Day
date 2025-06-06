#!/bin/bash
#SBATCH --job-name=xfer_I5d_eORCA12
#SBATCH --time=05:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_5d_3d_icemod_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 5d mean 
# icemod variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-10
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc

# Filepath to eORCA12 5d mean output files:
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

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Define the filepath to the 5d mean icemod files:
filepath_icemod=$filedir/1976/eORCA12_5d_icemod_1976*.nc
# Define the object prefix for the JASMIN object store:
object_prefix=I5d

echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 I5d variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":2160, "y":1803}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "Completed: Sent 1976 eORCA12 ERA-5 v1 I5d fields to JASMIN object store."

# -- Send eORCA12 5d mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1977..2024}
do
    # Define the filepath to the 5d mean icemod files:
    filepath_icemod=$filedir/${yr}/eORCA12_5d_icemod_${yr}*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 I5d variables to JASMIN object store..."
    msm_os update_with_dask -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                            -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                            -cs '{"x":2160, "y":1803}' \
                            -dc $dask_config_json -zv 3 || { echo "Error: msm_os update_with_dask failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 I5d fields to JASMIN object store."
done
