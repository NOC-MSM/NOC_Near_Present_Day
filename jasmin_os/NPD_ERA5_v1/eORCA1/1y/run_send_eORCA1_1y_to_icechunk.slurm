#!/bin/bash
#SBATCH --job-name=xfer_1y_eORCA1
#SBATCH --partition=compute
#SBATCH --time=04:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_1y_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA1 ERA-5 v1 annual-mean 
# output variables to Icechunk repos in the Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-05-22
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA1 ancillary file:
filepath_domain=/dssgfs01/scratch/npd/simulations/Domains/eORCA1/domain_cfg.nc

# Filepath to eORCA1 annual mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA1_ERA5_v1
filepath_gridT=$filedir/eORCA1_ERA5_1y_grid_T_*.nc
filepath_gridU=$filedir/eORCA1_ERA5_1y_grid_U_*.nc
filepath_gridV=$filedir/eORCA1_ERA5_1y_grid_V_*.nc
filepath_gridW=$filedir/eORCA1_ERA5_1y_grid_W_*.nc
filepath_icemod=$filedir/eORCA1_ERA5_1y_icemod_*.nc
filepath_scalar=$filedir/eORCA1_ERA5_1y_scalar_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=/dssgfs01/working/otooth/AtlantiS/jasmin_os/credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=/dssgfs01/scratch/otooth/NOC_NPD/20250512_eORCA1_ERA5v1/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA1/1y/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca1-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA1 ERA-5 v1 annual mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA1 ERA-5 v1 T1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":360,"y":331,"deptht":25}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA1 ERA5v1 T1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 U1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y \
                     -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                     -cs '{"x":360,"y":331,"depthu":25}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA1 ERA5v1 U1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 V1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y \
                     -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                     -cs '{"x":360,"y":331,"depthv":25}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA1 ERA5v1 V1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 W1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1y \
                      -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":360,"y":331,"depthw":25}' -dc $dask_config_json \
                      -br "main" -cm "Added eORCA1 ERA5v1 W1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 I1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":360,"y":331}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA1 ERA5v1 I1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 S1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y \
                     -dc $dask_config_json \
                     -br "main" -cm "Added eORCA1 ERA5v1 S1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA1 ERA-5 v1 annual mean fields to Icechunk repos."
