#!/bin/bash
#SBATCH --job-name=xfer_1y_eORCA025
#SBATCH --partition=compute
#SBATCH --time=05:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_1y_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA025 monthly-mean 
# output variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
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
filepath_icemod=$filedir/eORCA025_ERA5_1y_icemod_*.nc
filepath_scalar=$filedir/eORCA025_ERA5_1y_scalar_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Bucket and object prefix:
bucket=npd-eorca025-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 annual mean outputs to JASMIN OS -- #
echo "In Progress: Sending eORCA025 T1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1y \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":720, "y":603, "deptht":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/T1y/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/T1y/"}' \
                      -dlc '{"n_workers":40,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA025 U1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1y \
                      -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                      -cs '{"x":720, "y":603, "depthu":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/U1y/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/U1y/"}' \
                      -dlc '{"n_workers":40,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA025 V1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1y \
                      -gf $filepath_grid -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                      -cs '{"x":720, "y":603, "depthv":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/V1y/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/V1y/"}' \
                      -dlc '{"n_workers":40,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA025 W1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1y \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":720, "y":603, "depthw":25}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/W1y/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/W1y/"}' \
                      -dlc '{"n_workers":40,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA025 I1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_icemod -c $store_credentials_json -b $bucket -p I1y \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -cs '{"x":720, "y":603}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/I1y/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/I1y/"}' \
                      -dlc '{"n_workers":40,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "In Progress: Sending eORCA025 S1y variables to JASMIN object store..."
msm_os send_with_dask -f $filepath_scalar -c $store_credentials_json -b $bucket -p S1y \
                      -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/S1y/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/S1y/"}' \
                      -dlc '{"n_workers":10,"threads_per_worker":1,"memory_limit":"2GB"}'

echo "Completed: Sent eORCA025 annual mean fields to JASMIN object store."
