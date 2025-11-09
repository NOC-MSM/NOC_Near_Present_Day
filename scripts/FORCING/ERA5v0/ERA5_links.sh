#!/bin/bash

set -ex

linksdir='/dssgfs01/scratch/npd/forcing/ERA5/preprocessed/all_fields'
mkdir -p $linksdir
cd $linksdir

# Link all files
ln -s /dssgfs01/scratch/npd/forcing/ERA5/preprocessed/????/*/*nc .

# Link in weights files
ln -s /dssgfs01/working/atb299/NEMO_cfgs/ERA5_wgts/*nc .

# Fix date format
rename _19 _y19 *nc
rename _20 _y20 *nc
rename - m *nc

# Rename with variable names
rename 10m_u_component_of_wind u10 10m_u_component_of_wind*nc
rename 10m_v_component_of_wind v10 10m_v_component_of_wind*nc
rename 2m_dewpoint_temperature d2m 2m_dewpoint_temperature*nc
rename 2m_temperature t2m 2m_temperature*nc
rename mean_sea_level_pressure msl mean_sea_level_pressure*nc
rename mean_surface_downward_short_wave_radiation_flux msdwswrf mean_surface_downward_short_wave_radiation_flux*nc
rename mean_surface_downward_long_wave_radiation_flux msdwlwrf mean_surface_downward_long_wave_radiation_flux*nc
rename mean_total_precipitation_rate mtpr mean_total_precipitation_rate*nc
rename mean_snowfall_rate msr mean_snowfall_rate*nc

# Rename with variable names (wave data)
rename mean_wave_direction mwd mean_wave_direction*nc
rename mean_wave_period mwp mean_wave_period_y*nc
rename mean_wave_period_based_on_first_moment mp1 mean_wave_period_based_on_first_moment_y*nc
rename significant_height_of_combined_wind_waves_and_swell swh significant_height_of_combined_wind_waves_and_swell*nc
rename u_component_stokes_drift ust u_component_stokes_drift*nc
rename v_component_stokes_drift vst v_component_stokes_drift*nc

rename significant_height_of_total_swell shts significant_height_of_total_swell*nc
rename significant_height_of_wind_waves shww significant_height_of_wind_waves*nc
rename mean_wave_period_based_on_first_moment_for_swell p1ps mean_wave_period_based_on_first_moment_for_swell*nc
rename mean_wave_period_based_on_first_moment_for_wind_waves p1ww mean_wave_period_based_on_first_moment_for_wind_waves*nc
rename mean_direction_of_total_swell mdts mean_direction_of_total_swell*nc
rename mean_direction_of_wind_waves mdww mean_direction_of_wind_waves*nc
rename mean_period_of_total_swell mpts mean_period_of_total_swell*nc
rename mean_period_of_wind_waves mpww mean_period_of_wind_waves*nc

rename friction_velocity zust friction_velocity*nc
rename normalized_stress_into_ocean tauoc normalized_stress_into_ocean*nc

cd -
