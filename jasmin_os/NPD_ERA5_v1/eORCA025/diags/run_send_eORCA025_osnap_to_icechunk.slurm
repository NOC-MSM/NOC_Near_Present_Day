#!/bin/bash
#SBATCH --job-name=xfer_osnap_eORCA025
#SBATCH --partition=compute
#SBATCH --time=00:45:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_osnap_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA025 OSNAP section 
# variables for the RAPID, MOVE and SAMBA arrays Icechunk repos
# in the JASMIN Object Store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-08-15
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA025 OSNAP .nc file:
filedir=/dssgfs01/scratch/npd/diagnostics/eORCA025_ERA5_v1/osnap_ERA5_v1
filepath_osnap=$filedir/eORCA025_ERA5_v1_osnap_1976-01_2024-12.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=.../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send osnap eORCA025 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA025 ERA-5 v1 monthly OSNAP array outputs to Icechunk..."
ods send_to_icechunk -f $filepath_osnap -c $store_credentials_json -b $bucket -p "M1m/OSNAP" \
                     -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025 ERA5v1 OSNAP 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA025 ERA-5 v1 monthly OSNAP outputs to Icechunk repo."
