#!/bin/bash
#SBATCH --job-name=xfer_I5d_eORCA12
#SBATCH --partition=compute
#SBATCH --time=06:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA12_5d_icemod_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 5-day mean 
# SI3 ice model variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-03-03
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/work/n01/n01/atb299/NPD_eORCA12_v4.2/domain_cfg.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY

# Filepath to JASMIN OS credentials:
store_credentials_json=..../jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=..../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA12 5-day mean outputs to JASMIN OS -- #
for yr in {1976..2024}
do
    # Define the filepath to the 5-day mean I-point files:
    filepath_icemod=$filedir/eORCA12_ERA5_5d_icemod_${yr}*.nc
    # Define the object prefix for the JASMIN object store:
    object_prefix=I5d/${yr}

    echo "In Progress: Sending $yr eORCA12 ERA-5 v1 I5d variables to JASMIN object store..."
    msm_os send_with_dask -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                          -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                          -cs '{"x":720, "y":603}' \
                          -dc $dask_config_json -zv 3 || { echo "Error: msm_os update_with_dask failed."; exit 1; }

    echo "Completed: Sent $yr eORCA12 ERA-5 v1 5-day mean fields to JASMIN object store."
done
