#!/bin/bash
#SBATCH --job-name=xfer_S1m_eORCA12
#SBATCH --time=05:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_1d_S1m_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# 1-D scalar variables to Icechunk.
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

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the object prefix for the icechunk store:
object_prefix=S1m_1d

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over years:
for yr in {1980..2024}
do
    # Define the filepath to the monthly mean scalar files:
    filepath_scalar=$filedir/${yr}/eORCA12_1m_scalar_${yr}*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 S1m_1d variables to icechunk store..."
    ods update_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p $object_prefix \
                        -dc $dask_config_json || { echo "Error: ods update_icechunk failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 S1m_1d fields to icechunk store."
done
