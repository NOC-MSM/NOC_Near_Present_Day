#!/bin/bash
#SBATCH --job-name=xfer_1y_eORCA12
#SBATCH --time=05:00:00
#SBATCH --ntasks=30
#SBATCH --mem=75GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_1y_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 monthly-mean 
# output variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc

# Filepath to eORCA12 annual mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Define filepaths:
filepath_gridT=$filedir/????/eORCA12_1y_grid_T_*.nc
filepath_gridU=$filedir/????/eORCA12_1y_grid_U_*.nc
filepath_gridV=$filedir/????/eORCA12_1y_grid_V_*.nc
filepath_gridW=$filedir/????/eORCA12_1y_grid_W_*.nc
filepath_icemod=$filedir/????/eORCA12_1y_icemod_*.nc
filepath_scalar=$filedir/????/eORCA12_1y_scalar_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 annual mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA12 T1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":540, "y":721, "deptht":25}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA12 U1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                      -cs '{"x":540, "y":721, "depthu":25}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA12 V1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                      -cs '{"x":540, "y":721, "depthv":25}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA12 W1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":540, "y":721, "depthw":25}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA12 I1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":540, "y":721}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA12 S1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "Completed: Sent eORCA12 annual mean fields to JASMIN object store."
