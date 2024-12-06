#!/bin/bash

# -------------------------------------------------------------------------
# download_woa23_climatologies.sh
#
# Description: Download the World Ocean Atlas 2023 30-year climatologies
# using URLs provided by the National Centers for Environmental Information.
# WOA23 provides Climate Normals for temperature and salinity at 0.25 degree
# resolution as 30-year averages for initialising models.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2024-11-21
#
# -------------------------------------------------------------------------

echo "--- Downloading WOA23 climatologies ---"
echo "-> Downloading Upper 1500m WOA23 temperature climatologies"
# for i in {01..12}; do echo $i; wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/temperature/netcdf/decav71A0/0.25/woa23_decav71A0_t${i}_04.nc ; done

for i in {01..12}; do echo $i; wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/temperature/netcdf/decav81B0/0.25/woa23_decav81B0_t${i}_04.nc ; done

for i in {01..12}; do echo $i; wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/temperature/netcdf/decav91C0/0.25/woa23_decav91C0_t${i}_04.nc ; done
echo "Completed: Downloaded Upper 1500m WOA23 monthly temperature climatologies"

echo "-> Downloading Upper 1500m WOA23 monthly salinity climatologies"
for i in {01..12}; do echo $i; wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/salinity/netcdf/decav71A0/0.25/woa23_decav71A0_s${i}_04.nc ; done

for i in {01..12}; do echo $i; wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/salinity/netcdf/decav81B0/0.25/woa23_decav81B0_s${i}_04.nc ; done

for i in {01..12}; do echo $i; wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/salinity/netcdf/decav91C0/0.25/woa23_decav91C0_s${i}_04.nc ; done
echo "Completed: Downloaded Upper 1500m WOA23 monthly salinity climatologies"

# Winter average of full time period, covering full ocean depth:
echo "-> Downloading Full Depth WOA23 temperature and salinity climatologies..."
wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/temperature/netcdf/decav/0.25/woa23_decav_t13_04.nc
wget https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DATA/salinity/netcdf/decav/0.25/woa23_decav_s13_04.nc
echo "Completed: Downloaded Full Depth WOA23 temperature and salinity climatologies"
