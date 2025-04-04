#!/bin/bash
#SBATCH --job-name=xfer_domain_eORCA1
#SBATCH --partition=test
#SBATCH --time=00:15:00
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
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca1-era5v1

# Variables - save entire dataset to single zarr store - and chunks:
variables=consolidated

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send ancillary files to JASMIN OS -- #
echo "In Progress: Sending eORCA1 ERA-5 v1 mesh mask file to JASMIN object store..."
object_prefix=domain/domain_cfg
msm_os send_with_dask -f $filepath_mesh_mask -c $store_credentials_json -b $bucket -p $object_prefix \
                      -v $variables -cs '{"x":360,"y":331}' -zv 3 \
                      -dc $dask_config_json

echo "In Progress: Sending eORCA1 ERA-5 v1 sub-basins mask file to JASMIN object store..."
object_prefix=domain/subbasins
msm_os send_with_dask -f $filepath_subbasins -c $store_credentials_json -b $bucket -p $object_prefix \
                      -v $variables -cs '{"x":360,"y":331}' -zv 3 \
                      -dc $dask_config_json

echo "Completed: Sent eORCA1 ERA-5 v1 ancillary files to JASMIN object store."