# NOC NPD: Preparing ERA-5 Atmospheric Forcing.

This is a brief guide on how to prepare the bias corrected ERA-5 atmospheric forcing dataset used as the surface boundary condition in the NPD eORCA1, eORCA025 and eORCA12 simulations.

## Quick start on {Archer2|Anemone}

The following commands will check out and set up an instance of NPD. It is not advised to do this in your home directory.

```shell
git clone git@github.com:NOC-MSM/NOC_Near_Present_Day.git
cd NOC_Near_Present_Day
./setup {-s Archer2}
```

## Overview

To prepare ERA-5 atmospheric forcing fields, we follow the proceedure outlined below:

1. Download global ERA-5 atmospheric fields from the Copernicus Climate Data Store at hourly temporal resolution in monthly netCDF files.

2. Preprocess monthly netCDF files to perform flood filling over land, rechunking and variable renaming.

3. Perform bias correction of 2-m air temperature fields using ERA-5 sea ice concentration field and JRA-55 2-m air temperature climatology.

4. Create symbolic links to all processed ERA-5 atmospheric forcing files in the `/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/all_fields` directory.

In the following sections, we will discuss how to perform these steps in detail on the Anemone HPC.

## 1. Downloading ERA-5 Atmospheric Fields

Inside the `NOC_Near_Present_Day/scripts/FORCING/ERA5v1/` directory, users will find a bash script `run_download_ERA5_forcing.sh` to automate the downloading of ERA-5 data from the Copernicus Climate Data Store. 

One executing `run_download_ERA5_forcing.sh`, it calls the `download_ERA5_forcing.py` Python script, which (by default) downloads the following variables from the specified input `year` to present:

- 2m_temperature
- 2m_dewpoint_temperature
- 10m_u_component_of_wind
- 10m_v_component_of_wind
- mean_total_precipitation_rate
- mean_snowfall_rate
- mean_surface_downward_short_wave_radiation_flux
- mean_surface_downward_long_wave_radiation_flux
- mean_sea_level_pressure
- sea_ice_cover
- sea_surface_temperature

Each variable is downloaded for the entire globe at hourly temporal resolution & stored in a single netCDF file per month per variable.

These original ERA-5 variable netCDF files are stored according to their year in the `/dssgfs01/scratch/npd/forcing/ERA5/original/` directory and the `/dssgfs01/scratch/npd/forcing/ERA5/original_latest/` directory for the latest 3-months (since these are subject to change). 

## 2. Preprocessing ERA-5 Atmospheric Fields

Once we completed downloading our ERA-5 atmospheric forcing files using the `run_download_ERA5_forcing.sh` script, we next need to pre-process these fields.

There are three steps to pre-processing each ERA-5 variable netCDF file, all of which are perfomed by the the `run_preprocess_ERA5_forcing_original.slurm` and `run_preprocess_ERA5_forcing_original_latest.slurm` scripts.

1. All land grid points in the ERA-5 forcing field are flood filled using a land-sea mask created using the `create_land_sea_mask.py` script and the `ifthen` and `fillmiss3` operations available in the `cdo` library.

2. Atmospheric forcing netCDF files are rechunked using `ncks --cnk_dmn` to ensure longitude and latitude dimensions both have chunksizes of 24 values per chunk.

3. For atmospheric forcing files downloaded since 2025, the following variables are renamed to their previous names using `ncrename`

    - **mean_snowfall_rate**: avg_tsrwe >>> msr
    - **mean_surface_downward_short_wave_radiation_flux**:  avg_sdswrf >>> msdwswrf
    - **mean_surface_downward_long_wave_radiation_flux**:  avg_sdlwrf >>> msdwlwrf
    - **mean_total_precipitation_rate**: avg_tprate >>> mtpr

All pre-processed atmospheric forcing files are stored in the `/dssgfs01/scratch/npd/forcing/ERA5/preprocessed/` and `/dssgfs01/scratch/npd/forcing/ERA5/preprocessed_latest/` directories.

## 3. Perform Bias Correction of 2 m Air Temperature

There is a well established Surface Air Temperature (SAT) bias at high-latitudes in the ERA-5 atmospheric reanalysis (Tjernström & Graversen, 2009; Zampieri et al., 2023) owing to the poor representation of snow atop of sea ice (Batrak & Müller, 2019).

To account for these large biases, a climatological adjustment is applied to the ERA-5 hourly 2 m air temperature field over regions where (ERA-5) sea ice cover > 0%. 

Climatological offset factors are determined by calculating the difference between the long-term mean (1960-2019) monthly 2 m air temperature climatologies of ERA-5 and JRA55-do. These offsets are calculated using the `create_ERA5_t2m_si_adjustment.py` script available in the `ERA5v1` directory.

To avoid step-like transitions between monthly adjustments, 2 m air temperature offset factors are linearly interpolated in time to produce hourly surface forcing fields in the `apply_ERA5_t2m_si_adjustment.py` script.

To perform the 2 m air temperature bias correction for multiple years of preprocessed netCDF files, we can use the `run_apply_t2m_si_adjust_MY.slurm` script.

Alternatively, the `run_apply_t2m_si_adjust.slurm` script can be used when updating the ERA-5 atmospheric forcing dataset with less than 1 year of preprocessed files.

In both cases, the resulting monthly bias corrected 2 m air temperature netCDF files are stored in the `/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/` directory according to their year.

## 4. Create Links to ERA-5 Forcing Files

We first need to create symbolic links to all preprocessed ERA-5 atmospheric forcing files, excluding the bias corrected 2 m air temperature files.

**Note:** These commands should be performed inside the `/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/all_fields/` directory.

```bash
# Define NEMO forcing variables, excluding 2m_temperature
vars = ('2m_dewpoint_temperature' 'mean_sea_level_pressure' 'mean_surface_downward_long_wave_radiation_flux' 'mean_total_precipitation_rate' '10m_u_component_of_wind' '10m_v_component_of_wind' '2m_temperature' 'mean_snowfall_rate' 'mean_surface_downward_short_wave_radiation_flux')

# Create symbolic links:
for var in ${vars[*]}; do
    ln -s /dssgfs01/scratch/npd/forcing/ERA5/preprocessed/????/${var}/*nc .
done
```

For the bias corrected 2 m air temperature, we use the following path instead:

```bash
ln -s /dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/????/*.nc .
```

Next, we need to update the date format of the preprocessed links to be read by the surface boundary condition module in NEMO:

```bash
# Fix date format
rename _19 _y19 *nc
rename _20 _y20 *nc
rename - m *nc
```

Finally, we need to rename the variable names within the links to their shorter standard names:

```bash
# Rename with variable names
rename 10m_u_component_of_wind u10 10m_u_component_of_wind*nc
rename 10m_v_component_of_wind v10 10m_v_component_of_wind*nc
rename 2m_dewpoint_temperature d2m 2m_dewpoint_temperature*nc
rename 2m_temperature t2m 2m_temperature*nc
rename mean_sea_level_pressure msl mean_sea_level_pressure*nc
rename mean_snowfall_rate msr mean_snowfall_rate*nc
rename mean_surface_downward_long_wave_radiation_flux msdwlwrf mean_surface_downward_long_wave_radiation_flux*nc
rename mean_surface_downward_short_wave_radiation_flux msdwswrf mean_surface_downward_short_wave_radiation_flux*nc
rename mean_total_precipitation_rate mtpr mean_total_precipitation_rate*nc
```

We have now completed preparing the ERA-5 atmospheric forcing fields & we can perform a NOC Near-Present Day simulation by creating a further link to the `/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/all_fields` inside our NEMO run directory (this link will then be referenced directly in our `namelist_cfg` file).

## Contacts

For any issues regarding the preparation of ERA-5 atmospheric forcing data, contact Ollie Tooth (**oliver.tooth@noc.ac.uk**) or Adam Blaker (**atb299@noc.ac.uk**).
