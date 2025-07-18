#!/bin/bash
#SBATCH --job-name="xfer_eORCA12_icemod"
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --account=noc_msm
#SBATCH --partition=standard
#SBATCH --qos=standard

# ==============================================================
# run_send_eORCA12_1m_icemod_to_jasmin_os.slurm
#
# Description: SLURM script to send the eORCA12 monthly mean 
# icemod variables to the JASMIN object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-04-07
#
# ==============================================================
# -- Input arguments to msm-os -- #
# Filepath to eORCA12 ancillary file:
filepath_domain=/gws/nopw/j04/atlantis/otooth/NPD/eORCA12/domain_cfg.nc

# Filepath to eORCA12 monthly mean output files:
filedir=/gws/nopw/j04/atlantis/atb299/NPD/eORCA12

# Filepath to JASMIN OS credentials:
store_credentials_json=/home/users/o_tooth/NOC/jasmin_os/credentials/jasmin_object_store_credentials.json

# Bucket and object prefix:
bucket=npd-eorca12-era5v1

# -- Python Environment -- #
# Activate Python virtual environment:
module unload jaspy
source ~/miniforge3/bin/activate
conda activate env_jasmin_os

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Define the filepath to the monthly mean icemod files:
filepath_icemod=$filedir/1976/eORCA12_1m_icemod_1976*.nc
# Define the object prefix for the JASMIN object store:
object_prefix=I1m

echo "In Progress: Sending 1976 eORCA12 ERA-5 v1 I1m variables to JASMIN object store..."
msm_os send -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
            -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
            -cs '{"x":2160, "y":1803}' \
            -zv 3 || { echo "Error: msm_os send failed."; exit 1; }

echo "Completed: Sent 1976 eORCA12 ERA-5 v1 I1m fields to JASMIN object store."

# -- Send eORCA12 monthly mean outputs to JASMIN OS -- #
# Iterating over decades:
for yr in {1990..1996}
do
    # Define the filepath to the monthly mean icemod files:
    filepath_icemod=$filedir/${yr}/eORCA12_1m_icemod_${yr}*.nc

    echo "In Progress: Appending $yr eORCA12 ERA-5 v1 I1m variables to JASMIN object store..."
    msm_os update -f $filepath_icemod -c $store_credentials_json -b $bucket -p $object_prefix \
                  -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                  -cs '{"x":2160, "y":1803}' \
                  -zv 3 || { echo "Error: msm_os update_with_dask failed."; exit 1; }

    echo "Completed: Append $yr eORCA12 ERA-5 v1 I1m fields to JASMIN object store."
done
