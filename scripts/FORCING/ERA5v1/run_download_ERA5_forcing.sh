#!/bin/bash

# ------------------------------------------------------------
# run_download_ERA5_forcing.sh
# Description: Download ERA-5 atmospheric forcing data for a
# single year from the Copernicus Climate Data Store.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Last Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# ------------------------------------------------------------

# -- Load modules and activate conda environment -- #
module purge
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate /dssgfs01/working/otooth/conda_envs/env_npd_era5

# -- Input Parameters -- #
# Year to start ERA-5 monthly download from:
year=2025
# Directory to save historic ERA-5 monthly .nc files:
savedir=/dssgfs01/scratch/npd/forcing/ERA5/original/
# Directory to save latest 3-months of ERA-5 monthly .nc files:
latestdir=/dssgfs01/scratch/npd/forcing/ERA5/original_latest/
# Variable list to download - default is all variables:
# variables="2m_temperature 2m_dewpoint_temperature"

# -- Download ERA-5 Atmospheric Forcing -- #
echo "In Progress: Downloading ERA-5 atmospheric forcing data from ${year} to present..."

# Using default variable list:
python3 download_ERA5_forcing.py --year ${year} --savedir ${savedir} --latestdir ${latestdir} || { echo "Error: Downloading ERA-5 atmospheric forcing data failed..."; exit 1; }

echo "Completed: Downloaded ERA-5 atmospheric forcing data from ${year} to present."
