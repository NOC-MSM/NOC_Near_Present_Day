#!/bin/bash
#SBATCH --job-name=xfer_3d_I1m_eORCA025
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_3d_I1m_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA025 monthly mean 
# 3-D icemod variables to icechunk repositories.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-05-26
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filepath_domain=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 monthly mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1
filepath_icemod=$filedir/eORCA025_ERA5_1m_icemod_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=/dssgfs01/working/otooth/AtlantiS/jasmin_os/credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=/dssgfs01/scratch/otooth/NOC_NPD/20250512_eORCA1_ERA5v1/NOC_Near_Present_Day/jasmin_os/NPD_ERA5_v1/eORCA025/1m/dask_tmp/I1m/dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1
prefix=I1m

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 monthly mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA025 ERA-5 v1 I1m variables to Icechunk..."
ods send_to_icechunk -f $filepath_icemod -c $store_credentials_json -b $bucket -p $prefix \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":1440, "y":1206}' -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 I1m 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA025 ERA-5 v1 I1m monthly mean fields to Icechunk."
