#!/bin/bash
#SBATCH --job-name=xfer_metric_eORCA1
#SBATCH --partition=compute
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ==============================================================
# run_send_eORCA1_metric_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA1 METRIC section 
# variables for the RAPID, MOVE and SAMBA arrays to the JASMIN
# object store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2025-08-15
#
# ==============================================================
# -- Input arguments to ods -- #
# Filepaths to eORCA1 ERA-5 v1 diagnostics:
filedir=/dssgfs01/scratch/npd/diagnostics/eORCA1_ERA5_v1/metric_ERA5_v1
filepath_rapid=$filedir/eORCA1_ERA5_v1_1976-01-2024-12_natl_meridional_transports_at_26N.nc
filepath_move=$filedir/eORCA1_ERA5_v1_1976-01-2024-12_natl_meridional_transports_at_16N.nc
filepath_samba=$filedir/eORCA1_ERA5_v1_1976-01-2024-12_natl_meridional_transports_at_34_5S.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=.../credentials/jasmin_os_credentials.json

# Filepath to Dask Configuration file:
dask_config_json=.../dask_config.json

# Bucket and object prefix:
bucket=npd-eorca1-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# -- Send METRIC eORCA1 ERA-5 v1 monthly mean diagnostics to JASMIN OS -- #
echo "In Progress: Sending eORCA1 ERA-5 v1 monthly METRIC RAPID 26.5N outputs to Icechunk..."
ods send_to_icechunk -f $filepath_rapid -c $store_credentials_json -b $bucket -p "M1m/RAPID_26N" \
                     -dc $dask_config_json -ad "time" \
                     -br "main" -cm "Added eORCA1 ERA5v1 METRIC RAPID 26N 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 monthly METRIC MOVE 16N outputs to Icechunk..."
ods send_to_icechunk -f $filepath_move -c $store_credentials_json -b $bucket -p "M1m/MOVE_16N" \
                     -dc $dask_config_json -ad "time" \
                     -br "main" -cm "Added eORCA1 ERA5v1 METRIC MOVE 16N 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "In Progress: Sending eORCA1 ERA-5 v1 monthly METRIC SAMBA 34.5S outputs to Icechunk..."
ods send_to_icechunk -f $filepath_samba -c $store_credentials_json -b $bucket -p "M1m/SAMBA_34_5S" \
                     -dc $dask_config_json -ad "time" \
                     -br "main" -cm "Added eORCA1 ERA5v1 SAMBA 34.5S 1976-2024" || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA1 ERA-5 v1 monthly METRIC outputs to Icechunk repos."
