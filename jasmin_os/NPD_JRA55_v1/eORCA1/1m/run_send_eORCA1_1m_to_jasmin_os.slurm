#!/bin/bash
#SBATCH --job-name=xfer_1m_eORCA1
#SBATCH --partition=compute
#SBATCH --time=05:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_1m_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA1 monthly-mean 
# output variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-01-02
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA1 ancillary file:
filepath_grid=/dssgfs01/scratch/otooth/npd_data/simulations/Domains/eORCA1/domain_cfg.nc

# Filepath to eORCA1 monthly mean output files:
filedir=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA1_JRA55/exp_npd_v1
filepath_gridT=$filedir/eORCA1_1m_grid_T_*.nc
filepath_gridU=$filedir/eORCA1_1m_grid_U_*.nc
filepath_gridV=$filedir/eORCA1_1m_grid_V_*.nc
filepath_gridW=$filedir/eORCA1_1m_grid_W_*.nc
filepath_scalar=$filedir/eORCA1_1m_scalar_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Bucket and object prefix:
bucket=npd-eorca1-jra55v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA1 monthly mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA1 T1m variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1m \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":360,"y":331,"deptht":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/"}' \
                      -dlc '{"n_workers":30,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA1 U1m variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1m \
                      -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                      -cs '{"x":360,"y":331,"depthu":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/"}' \
                      -dlc '{"n_workers":30,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA1 V1m variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1m \
                      -gf $filepath_grid -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                      -cs '{"x":360,"y":331,"depthv":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/"}' \
                      -dlc '{"n_workers":30,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA1 W1m variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1m \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":360,"y":331,"depthw":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/"}' \
                      -dlc '{"n_workers":30,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA1 scalar variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1m \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/"}' \
                      -dlc '{"n_workers":5,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "Completed: Sent eORCA1 monthly mean fields to JASMIN object store."
