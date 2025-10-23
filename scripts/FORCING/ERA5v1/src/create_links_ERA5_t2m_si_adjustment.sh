#!/bin/bash

set -ex
# ------------------------------------------------------------
# create_links_ERA5_t2m_si_adjustment.sh
# Description: Creates links to adjusted ERA5 2m temperature
# files in a directory.
# ------------------------------------------------------------
linksdir='/dssgfs01/scratch/npd/forcing/ERA5/preprocessed/all_fields'
mkdir -p $linksdir
cd $linksdir

# -- Link all ERA-5 forcing files -- #
ln -s /dssgfs01/scratch/npd/forcing/ERA5/preprocessed/????/*/*nc .

# -- Link all ERA-5 weight files -- #
ln -s /dssgfs01/working/atb299/NEMO_cfgs/ERA5_wgts/*nc .

# -- Update date format -- #
rename _19 _y19 *nc
rename _20 _y20 *nc
rename - m *nc

# -- Rename with variable names -- #
rename 10m_u_component_of_wind u10 10m_u_component_of_wind*nc
rename 10m_v_component_of_wind v10 10m_v_component_of_wind*nc
rename 2m_dewpoint_temperature d2m 2m_dewpoint_temperature*nc
rename 2m_temperature t2m 2m_temperature*nc
rename mean_sea_level_pressure msl mean_sea_level_pressure*nc
rename mean_surface_downward_short_wave_radiation_flux msdwswrf mean_surface_downward_short_wave_radiation_flux*nc
rename mean_surface_downward_long_wave_radiation_flux msdwlwrf mean_surface_downward_long_wave_radiation_flux*nc
rename mean_total_precipitation_rate mtpr mean_total_precipitation_rate*nc
rename mean_snowfall_rate msr mean_snowfall_rate*nc

cd -