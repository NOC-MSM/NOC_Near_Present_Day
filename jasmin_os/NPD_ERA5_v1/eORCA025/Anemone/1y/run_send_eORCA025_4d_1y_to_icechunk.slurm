#!/bin/bash
#SBATCH --job-name=xfer_4d_1y_eORCA025
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_4d_1y_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA025 annual-mean 
# 4-D output variables to the Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-05-26
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 annual mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1
filepath_gridT=$filedir/eORCA025_ERA5_1y_grid_T_*.nc
filepath_gridU=$filedir/eORCA025_ERA5_1y_grid_U_*.nc
filepath_gridV=$filedir/eORCA025_ERA5_1y_grid_V_*.nc
filepath_gridW=$filedir/eORCA025_ERA5_1y_grid_W_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=/dssgfs01/working/otooth/AtlantiS/jasmin_os/credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=/dssgfs01/scratch/otooth/NOC_NPD/20250512_eORCA1_ERA5v1/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA025/1y/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 annual mean outputs to JASMIN OS -- #
variables_gridT="time_centered_bounds time_counter_bounds e3t thetao_con so_abs rsdo vohflisf vohfcisf vowflisf ocontempadvect ocontemppmdiff ocontempdiff ttrd_evd_li rsdoabsorb ttrd_bbl_li ttrd_atf_li ocontemptend osaltadvect osaltpmdiff osaltdiff strd_evd_li osalttend strd_bbl_li strd_atf_li"

echo "In Progress: Sending eORCA025 4-D T1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y_4d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":720, "y":603, "deptht":25}' -v $variables_gridT  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 T1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }
                    

variables_gridU="time_centered_bounds time_counter_bounds e3u uo uo_eiv u2o umo"
echo "In Progress: Sending eORCA025 4-D U1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y_4d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                     -cs '{"x":720, "y":603, "depthu":25}' -v $variables_gridU  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 U1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

variables_gridV="time_centered_bounds time_counter_bounds e3v vo vo_eiv v2o vmo"
echo "In Progress: Sending eORCA025 4-D V1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y_4d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                     -cs '{"x":720, "y":603, "depthv":25}' -v $variables_gridV  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 V1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA025 4-D W1y variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1y_4d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":720, "y":603, "depthw":25}'  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 W1y 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA025 ERA-5 v1 4-D annual mean fields to Icechunk."
