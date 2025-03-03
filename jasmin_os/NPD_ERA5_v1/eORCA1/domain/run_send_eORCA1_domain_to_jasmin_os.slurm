#!/bin/bash
#SBATCH --job-name=xfer_domain_eORCA1
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_domain_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA1 model domain 
# variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA1 ancillary files:
filepath_mesh_mask=/dssgfs01/scratch/npd/simulations/Domains/eORCA1/mesh_mask.nc
filepath_subbasins=/dssgfs01/scratch/npd/simulations/Domains/eORCA1/subbasins_CMIP6.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Bucket and object prefix:
bucket=npd-eorca1-era5v1
object_prefix=domain

# Variables - save entire dataset to single zarr store - and chunks:
variables=compact

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send ancillary files to JASMIN OS -- #
echo "In Progress: Sending eORCA1 ERA-5 v1 mesh mask file to JASMIN object store..."
msm_os send_with_dask -f $filepath_mesh_mask -c $store_credentials_json -b $bucket -p $object_prefix -v $variables \
                      -cs '{"x":360,"y":331}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/"}' \
                      -dlc '{"n_workers":10,"threads_per_worker":4,"memory_limit":"3GB"}'

echo "In Progress: Sending eORCA1 ERA-5 v1 sub-basins mask file to JASMIN object store..."
msm_os send_with_dask -f $filepath_subbasins -c $store_credentials_json -b $bucket -p $object_prefix -v $variables \
                      -cs '{"x":360,"y":331}' \
                      -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/jasmin_os_tmp/"}' \
                      -dlc '{"n_workers":10,"threads_per_worker":4,"memory_limit":"3GB"}'

echo "Completed: Sent eORCA1 ERA-5 v1 ancillary files to JASMIN object store."