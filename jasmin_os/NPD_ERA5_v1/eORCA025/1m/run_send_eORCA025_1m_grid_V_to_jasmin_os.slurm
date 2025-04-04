#!/bin/bash
#SBATCH --job-name=xfer_V_eORCA025
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_1m_grid_V_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA025 monthly mean 
# V-point variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filepath_domain=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 monthly mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1
filepath_gridV=$filedir/eORCA025_ERA5_1m_grid_V_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 monthly mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA025 ERA-5 v1 V1m variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1m \
                      -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                      -cs '{"x":720, "y":603, "depthv":25}' -zv 3 \
                      -dc $dask_config_json || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "Completed: Sent eORCA025 ERA-5 v1 monthly mean fields to JASMIN object store."
