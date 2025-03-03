#!/bin/bash
#SBATCH --job-name=xfer_T_eORCA025
#SBATCH --partition=compute
#SBATCH --time=10:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_1m_grid_T_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA025 monthly mean 
# T-point variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 monthly mean output files:
filedir=/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1
filepath_gridT=$filedir/eORCA025_ERA5_1m_grid_T_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 monthly mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA025 ERA-5 v1 T1m variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1m \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":720, "y":603, "deptht":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/T1m/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/T1m/"}' \
                      -dlc '{"n_workers":40,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "Completed: Sent eORCA025 ERA-5 v1 monthly mean fields to JASMIN object store."
