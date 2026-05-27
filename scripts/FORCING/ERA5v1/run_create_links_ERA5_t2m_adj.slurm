#!/bin/bash
#SBATCH --job-name=create_era5_t2m_adj_links
#SBATCH --time=00:05:00
#SBATCH --partition=test
#SBATCH --ntasks=1

# ------------------------------------------------------------
# run_create_links_ERA5_t2m_adj.slurm
# Description: Create symbolic links for climatologically
# adjusted ERA-5 2m temperature files.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Last Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# ------------------------------------------------------------
set -euo pipefail

# ==== Input Parameters ==== #
# Year to prepare symbolic links for:
year=2026

# Input directory containing preprocessed ERA-5 monthly .nc files:
inputdir="/dssgfs01/scratch/npd/forcing/ERA5/preprocessed_latest"

# Directory to save historic ERA-5 monthly .nc files:
linksdir="/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/all_fields/"

# === Create Symbolic Links for Climatologically Adjusted ERA-5 2m Temperature === #
echo "In Progress: Creating symbolic links for climatologically adjusted ERA-5 2m temperature files for ${year} ..."

# Move to directory of symbolic links for ERA-5 adjusted forcing:
mkdir -p $linksdir
cd $linksdir

# Define NEMO forcing variables, excluding 2m_temperature:
vars="2m_dewpoint_temperature mean_sea_level_pressure mean_surface_downward_long_wave_radiation_flux mean_total_precipitation_rate 10m_u_component_of_wind 10m_v_component_of_wind mean_snowfall_rate mean_surface_downward_short_wave_radiation_flux"

# Create symbolic links for preprocessed ERA-5 variables:
for var in ${vars}; do
    ln -s ${inputdir}/${year}/${var}/*nc .
done

# Create symbolic links for climatologically adjusted ERA-5 2m temperature:
ln -s /dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/${year}/2m_temperature*.nc .

# === Update File Naming Format === #
# Update ERA-5 date format (i.e., _2025-01 -> _y2025m01):
# rename _19 _y19 *nc
rename _20 _y20 *nc
rename - m *nc

# Update ERA-5 variable names to short format:
rename 10m_u_component_of_wind u10 10m_u_component_of_wind*nc
rename 10m_v_component_of_wind v10 10m_v_component_of_wind*nc
rename 2m_dewpoint_temperature d2m 2m_dewpoint_temperature*nc
rename 2m_temperature t2m 2m_temperature*nc
rename mean_sea_level_pressure msl mean_sea_level_pressure*nc
rename mean_snowfall_rate msr mean_snowfall_rate*nc
rename mean_surface_downward_long_wave_radiation_flux msdwlwrf mean_surface_downward_long_wave_radiation_flux*nc
rename mean_surface_downward_short_wave_radiation_flux msdwswrf mean_surface_downward_short_wave_radiation_flux*nc
rename mean_total_precipitation_rate mtpr mean_total_precipitation_rate*nc

echo "Completed: Created symbolic links for climatologically adjusted ERA-5 2m temperature files for ${year}."
cd -
