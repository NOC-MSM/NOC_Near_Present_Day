#!/bin/bash
#SBATCH --job-name=xfer_osnap_eORCA1
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_osnap_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA1 OSNAP section 
# variables for the RAPID, MOVE and SAMBA arrays Icechunk repos
# in the JASMIN Object Store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-05-2
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepath to eORCA1 OSNAP .nc file:
filedir=/dssgfs01/scratch/npd/simulations/eORCA1_ERA5_v1/osnap_ERA5_v1
filepath_osnap=$filedir/eORCA1_ERA5_1m_osnap_1976_2024.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca1-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send osnap eORCA1 monthly mean output to JASMIN OS -- #
echo "In Progress: Sending eORCA1 ERA-5 v1 monthly OSNAP array outputs to Icechunk..."
ods send_to_icechunk -f $filepath_osnap -c $store_credentials_json -b $bucket -p M1m/OSNAP \
                     -dc $dask_config_json \
                     -br "main" -cm "Added eORCA1 ERA5v1 OSNAP 1976-2024" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA1 ERA-5 v1 monthly OSNAP outputs to Icechunk repo."
