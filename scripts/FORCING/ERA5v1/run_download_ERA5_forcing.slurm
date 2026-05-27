#!/bin/bash
#SBATCH --job-name=download_era5_forcing
#SBATCH --time=06:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --exclusive

# ------------------------------------------------------------
# run_download_ERA5_forcing.sh
# Description: Download ERA-5 atmospheric forcing data for a
# single year from the Copernicus Climate Data Store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Last Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# ------------------------------------------------------------
set -euo pipefail

# -- Load modules and activate conda environment -- #
module purge
source /dssgfs01/working/otooth/miniforge3/bin/activate
conda activate /dssgfs01/working/otooth/conda_envs/env_npd_era5

# ==== Input Parameters ==== #
# Year to start ERA-5 monthly download from:
year=2026

# Optional: End date to include partial month download:
end_date="2026-06-01"

# Directory to save historic ERA-5 monthly .nc files:
savedir=/dssgfs01/scratch/npd/forcing/ERA5/original/

# Directory to save latest 3-months of ERA-5 monthly .nc files:
latestdir=/dssgfs01/scratch/npd/forcing/ERA5/original_latest/

# Variable list to download - default is all variables:
# variables="2m_temperature 2m_dewpoint_temperature"

# -- Download ERA-5 Atmospheric Forcing -- #
echo "In Progress: Downloading ERA-5 atmospheric forcing data from ${year} to present..."

# Using default variable list:
# python3 download_ERA5_forcing.py --year ${year} --savedir ${savedir} --latestdir ${latestdir}
python3 download_ERA5_forcing.py --year ${year} --savedir ${savedir} --latestdir ${latestdir} --enddate ${end_date}

echo "Completed: Downloaded ERA-5 atmospheric forcing data from ${year} to present."
