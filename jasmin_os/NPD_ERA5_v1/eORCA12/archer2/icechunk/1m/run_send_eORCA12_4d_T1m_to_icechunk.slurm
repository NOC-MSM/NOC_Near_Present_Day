#!/bin/bash
#SBATCH --job-name=xfer_T1m_4d_eORCA12
#SBATCH --time=18:00:00
#SBATCH --ntasks=8
#SBATCH --mem=24GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_4d_T1m_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# T-point 4D variables to Icechunk.
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

# Define 4-dimensional variables:
variables_gridT="time_centered_bounds time_counter_bounds e3t thetao_con so_abs rsdo vohflisf vohfcisf vowflisf"

# Filepath to JASMIN OS credentials:
store_credentials_json=/work/n01/n01/otooth/NOC-NPD/credentials/jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA12/archer2/icechunk/1y/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the object prefix for the icechunk store:
object_prefix=T1m_4d

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over years:
for yr in {1980..2023}
do
    # Define the filepath to the monthly mean T-point files:
    filepath_gridT=$filedir/????/eORCA12_1m_grid_T_${yr}*.nc

    echo "In Progress: Sending ${yr} eORCA12 ERA-5 v1 T1m 4D variables to icechunk store..."
    ods update_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p $object_prefix \
                        -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                        -cs '{"x":540, "y":721, "deptht":25}' -v $variables_gridT -dc $dask_config_json \
                        -br "main" -cm "Added eORCA12 ERA5v1 T1m 1976-2024" || { echo "Error: ods update_icechunk failed."; exit 1; }

    echo "Completed: Sent ${yr} eORCA12 ERA-5 v1 T1m 4D variables to icechunk store."
done
