#!/bin/bash
#SBATCH --job-name=xfer_W5d_eORCA025
#SBATCH --partition=compute
#SBATCH --time=08:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_W5d_to_zarr.slurm
#
# Description: SLURM script to send the eORCA025 5-day mean 
# W-point variables to the Zarr store.
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

# -- Send eORCA025 5-day mean outputs to JASMIN OS -- #
for yr in {1991..2024}
do
    # Define the filepath to the 5-day mean W-point files:
    filepath_gridW=$filedir/eORCA025_ERA5_5d_grid_W_${yr}*.nc
    # Define the object prefix for the Zarr store:
    object_prefix=W5d/${yr}

    echo "In Progress: Sending $yr eORCA025 ERA-5 v1 W5d variables to Zarr store..."
    msm_os send_with_dask -f $filepath_gridW -c $store_credentials_json -b $bucket -p $object_prefix \
                          -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                          -cs '{"x":720, "y":603, "depthw":25}' \
                          -dc $dask_config_json -zv 3 || { echo "Error: msm_os send_with_dask failed."; exit 1; }

    echo "Completed: Sent $yr eORCA025 ERA-5 v1 5-day mean fields to Zarr store."
done
