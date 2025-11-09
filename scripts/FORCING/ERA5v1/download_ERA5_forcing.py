"""
download_era5_forcing.py

Description: This script downloads ERA-5 atmospheric forcing data
from the Copernicus Climate Data Store.

Created By: Jenny Mecking (jmecki@noc.ac.uk)
Last Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# --- Import Dependencies --- #
import os
import cdsapi
import logging
import datetime
import argparse
import numpy as np
from datetime import date, timedelta
from tenacity import retry, stop_after_attempt, wait_fixed

# -- Configure Argument Parser -- #
# Define the argument parser:
parser = argparse.ArgumentParser(description='Download ERA-5 atmospheric forcing data from the Copernicus Climate Data Store at 3-hour frequency in monthly netCDF output files.')
# Define input arguments:
parser.add_argument('-y','--year', help='Year to start ERA-5 monthly download from.', required=True, type=int, default=2025)
parser.add_argument('-s','--savedir', help='Directory to save historic ERA-5 monthly netCDF files.', required=True)
parser.add_argument('-l','--latestdir', help='Directory to save latest 3-months of ERA-5 monthly netCDF files.', required=True)
parser.add_argument('-v','--variables', help='List of ERA-5 variables to download.', required=False, nargs='+', default=['2m_temperature', '2m_dewpoint_temperature', '10m_u_component_of_wind', '10m_v_component_of_wind',
 'mean_total_precipitation_rate', 'mean_snowfall_rate', 'mean_surface_downward_short_wave_radiation_flux', 'mean_surface_downward_long_wave_radiation_flux', 'mean_sea_level_pressure', 'sea_ice_cover', 'sea_surface_temperature'])

# --- Configure Logging --- #
logging.basicConfig(
    filename="download_era5_forcing.log",
    encoding="utf-8",
    filemode="w",
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
    )

# -- Define User Inputs -- #
# Parse input arguments as dict:
args = vars(parser.parse_args())
# Store user inputs:
ini_year = args['year']
save_dir = args['savedir']
latest_dir = args['latestdir']
variables = args['variables']

# -- Define download_ERA5_data() -- #
@retry(stop=stop_after_attempt(3), wait=wait_fixed(10))
def download_ERA5_data(client, date:datetime.date, variable:str, outfile:str) -> None:
    """
    Download ERA-5 3-hourly data on single levels using the Copernicus
    Climate Data Store (CDS) API.

    Parameters:
    -----------
    client: cdsapi.Client
        Climate Data Store API client instance.
    date: date
        Date of the data to download.
    variable: str
        Variable to download.
    outfile: str
        Path to save the downloaded file.
    
    Returns:
    --------
    None
    """
    # -- Verify Inputs -- #
    if not isinstance(date, datetime.date):
        raise TypeError("date must be a date object.")
    if not isinstance(variable, str):
        raise TypeError("variable must be a string.")
    if not isinstance(outfile, str):
        raise TypeError("outfile must be a string.")
    
    # -- Download Data -- #
    # Download ERA-5 data as a NetCDF file:
    client.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': variable,
            'year': str(date.year),
            'month': "{month:02d}".format(month=date.month),
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'format': 'netcdf',
        },
        outfile)
    
    return

# -- Prepare Dates -- #
ini_date = date(year=ini_year, month=1, day=1)
final_date = date.today()
# Excluding current month - downloads latest full month of data:
dates_monthly = np.arange(ini_date, final_date, dtype='datetime64[M]').astype(datetime.date)

# -- Initialise CDS API -- #
client = cdsapi.Client()
# Define number of download retries:
n_retry = 3

# -- Download ERA5 Data -- #
for date_mo in dates_monthly:
    for variable in variables:
        # Define output directory according to date:
        if date_mo > (final_date - timedelta(weeks=12)):
            outdir = f"{latest_dir}{date_mo.year}/{variable}/"
        else:
            outdir = f"{save_dir}{date_mo.year}/{variable}/"

        # Define output filepath:
        outfile = (outdir + "{variable}_{year}-{month:02d}.nc".format(variable=variable, year=date_mo.year, month=date_mo.month))

        # Check if filepath exists:
        if os.path.isfile(outfile):
            logging.info(f"Skipping file: {outfile} already exists.")
        else:
            # Check if directory exists:
            if not os.path.exists(outdir):
                # Create a new file directory:
                os.makedirs(outdir)

            # Download ERA5 data with a maximum of 3 retries:
            logging.info(f"Downloading {variable} data for {date_mo.year}-{date_mo.month}.")
            download_ERA5_data(client, date_mo, variable, outfile)
