#!/bin/bash
#SBATCH --job-name=xfer_1d_eORCA025
#SBATCH --partition=compute
#SBATCH --time=03:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ===================================================================
# run_send_eORCA025_1d_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA025 ERA-5 v1 daily-mean 
# output variables to Icechunk repos in the Icechunk.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-05-22
#
# ===================================================================
# -- Input arguments to ods -- #
# Filepath to eORCA025 ancillary file:
filepath_domain=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 daily mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1
filepath_gridT=$filedir/eORCA025_ERA5_1d_grid_T_*.nc
filepath_gridU=$filedir/eORCA025_ERA5_1d_grid_U_*.nc
filepath_gridV=$filedir/eORCA025_ERA5_1d_grid_V_*.nc
filepath_scalar=$filedir/eORCA025_ERA5_1d_scalar_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=/dssgfs01/working/otooth/AtlantiS/jasmin_os/credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=/dssgfs01/scratch/otooth/NOC_NPD/20250512_eORCA025_ERA5v1/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA025/1d/dask_config.json

# Bucket and object prefix:
bucket=npd-eORCA025-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 ERA-5 v1 daily mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA025 ERA-5 v1 T1d variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":1440, "y":1206}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 T1d 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA025 ERA-5 v1 U1d variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                     -cs '{"x":1440, "y":1206}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 U1d 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA025 ERA-5 v1 V1d variables to Icechunk..."
ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                     -cs '{"x":1440, "y":1206}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 V1d 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA025 ERA-5 v1 scalar variables to Icechunk..."
ods send_to_icechunk -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1d \
                     -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 S1d 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA025 ERA-5 v1 daily mean fields to Icechunk repos."
