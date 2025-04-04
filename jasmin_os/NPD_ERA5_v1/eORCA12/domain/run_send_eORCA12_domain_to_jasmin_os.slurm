#!/bin/bash
#SBATCH --job-name=xfer_domain_eORCA12
#SBATCH --time=00:30:00
#SBATCH --ntasks=12
#SBATCH --mem=30GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA12_domain_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 model domain 
# variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-12
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_mesh_mask=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc
# filepath_subbasins=/dssgfs01/scratch/npd/simulations/Domains/subbasins.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# Variables - save entire dataset to single zarr store - and chunks:
variables=consolidated

# -- Python Environment -- #
# Activate Python virtual environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send ancillary files to JASMIN OS -- #
echo "In Progress: Sending eORCA12 ERA-5 v1 mesh mask file to JASMIN object store..."
object_prefix=domain/mesh_mask
msm_os send_with_dask -f $filepath_mesh_mask -c $store_credentials_json -b $bucket \
                      -p $object_prefix -v $variables \
                      -cs '{"x":540,"y":721,"nav_lev":25}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "In Progress: Sending eORCA12 ERA-5 v1 mesh mask file to JASMIN object store..."
object_prefix=domain/subbasins
msm_os send_with_dask -f $filepath_subbasins -c $store_credentials_json -b $bucket \
                      -p $object_prefix -v $variables \
                      -cs '{"x":540,"y":721}' \
                      -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

echo "Completed: Sent eORCA12 ERA-5 v1 ancillary files to JASMIN object store."