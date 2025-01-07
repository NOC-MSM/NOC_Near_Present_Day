#!/bin/bash
#SBATCH --job-name=xfer_U5d_eORCA025
#SBATCH --partition=compute
#SBATCH --time=06:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA025_5d_grid_U_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA025 5-day mean 
# U-point variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-01-02
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/dssgfs01/scratch/otooth/npd_data/simulations/Domains/eORCA025/domain_cfg.nc

# Filepath to eORCA025 monthly mean output files:
filedir=/dssgfs01/scratch/otooth/npd_data/simulations/eORCA025_JRA55/exp_npd_v1

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Bucket and object prefix:
bucket=npd-eorca025-jra55v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send eORCA025 5-day mean outputs to JASMIN OS -- #
for yr in {1976..2024}
do
    # Define the filepath to the 5-day mean U-point files:
    filepath_gridU=$filedir/eORCA025_5d_grid_U_${yr}*.nc
    # Define the object prefix for the JASMIN object store:
    object_prefix=U5d/${yr}

    echo "In Progress: Sending $yr eORCA025 U5d variables to JASMIN object store..."
    msm_os send_with_dask -f $filepath_gridU -c $store_credentials_json -b $bucket -p $object_prefix \
                          -gf $filepath_grid -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                          -cs '{"x":720, "y":603, "depthu":25}' \
                          -dco '{"temporary_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/gridU/","local_directory":"/dssgfs01/scratch/otooth/NOC_NPD/20250102/NOC_Near_Present_Day/jasmin_os/tmp/gridU/"}' \
                          -dlc '{"n_workers":35,"threads_per_worker":1,"memory_limit":"2.5GB"}'

    echo "Completed: Sent $yr eORCA025 5-day mean fields to JASMIN object store."
done
