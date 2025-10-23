"""
apply_ERA5_t2m_si_adjustment.py

Description: This script applies monthly climatological correction to ERA-5
3-hourly 2m air temperature where sea ice cover exists.

Created By: Adam Blaker (atb299@noc.ac.uk)
Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import Dependencies -- #
import os
import glob
import logging
import argparse
import xarray as xr

# -- Configure Argument Parser -- #
# Define the argument parser:
parser = argparse.ArgumentParser(description='Apply monthly climatological correction to ERA-5 2m air temperature at 3-hourly frequency above sea ice.')
# Define input arguments:
parser.add_argument('-y','--year', help='Year to apply monthly climatological correction.', type=int, required=True)
parser.add_argument('-o','--outdir', help='Directory to save climatologically corrected ERA-5 2m temperature netCDF files.', default="/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj/")

# --- Configure Logging --- #
logging.basicConfig(
    filename="apply_t2m_si_adjustment.log",
    encoding="utf-8",
    filemode="w",
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
    )

# -- Define format_ERA5_data() -- #
def get_encoding(ds: xr.Dataset, yr: int) -> dict:
    """
    Define encoding for output netCDF file based on input 2m
    air temperature dataset.

    Parameters:
    -----------
    ds: xr.Dataset
        Input 2m temperature dataset to define encoding for.
        Must contain variable ``t2m`` storing 2m temperature data.
    yr: int
        Year of the input dataset.
    
    Returns:
    --------
    encoding: dict
        Encoding for output netCDF file.
    """
    # -- Verify Inputs -- #
    if not isinstance(ds, xr.Dataset):
        raise ValueError("Input dataset must be an xarray Dataset.")
    if not isinstance(yr, int):
        raise ValueError("Year must be an integer.")
    if 't2m' not in ds.data_vars:
        raise ValueError("Input dataset must contain variable ``t2m`` storing 2m temperature data.")
    
    # -- Define Encoding -- #
    # Determine/set chunking and scale/offset.
    # Note ERA5 file formats changed in 2024.
    if yr < 2024:
        encoding={'t2m': {'chunksizes': (1, 24, 24), 
                          'dtype': 'int16', 
                          'zlib': True, 
                          'complevel': 1,
                          'missing_value': -32767,
                          '_FillValue': -32767,
                          'scale_factor': ds['t2m'].encoding['scale_factor'],
                          'add_offset': ds['t2m'].encoding['add_offset']
                          }}
    else:
        encoding={'t2m': {'chunksizes': (1, 24, 24),
                          'dtype': 'float32',
                          'zlib': True,
                          'complevel': 1}}

    return encoding

# -- Define User Inputs -- #
# Parse input arguments as dict:
args = vars(parser.parse_args())
# Define year to apply monthly climatological correction:
yr = args["year"]
# Define output directory:
outdir = args["outdir"]
if outdir.endswith("/"):
    outdir = f"{outdir}{yr}"
else:
    outdir = f"{outdir}/{yr}"

# -- Define File Paths -- #
# ERA-5 2m temperature input directory:
t2m_fpath = f"/dssgfs01/scratch/npd/forcing/ERA5/preprocessed/{yr}/2m_temperature/*.nc"
t2m_files = sorted(glob.glob(t2m_fpath))
n_files = len(t2m_files)

# ERA-5 sea ice cover input directory:
sic_fpath = f"/dssgfs01/scratch/npd/forcing/ERA5/original/{yr}/sea_ice_cover/*.nc"
sic_files = sorted(glob.glob(sic_fpath))

# Check number of t2m and sic files match:
if (n_files != len(sic_files)):
    raise ValueError("Unequal number of 2m temperature and sea ice cover netCDF files.")

# Check if output directory exists:
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Define month names:
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# -- Apply Climatological Correction -- #
# Iterate over monthly 2m temperature and sea ice cover files:
for n in range(n_files):
    # Open 2m temperature dataset:   
    ds_t2m = xr.open_dataset(t2m_files[n])
    if ("valid_time" in ds_t2m.coords):
        ds_t2m = ds_t2m.rename({"valid_time": "time"})
    logging.info(f"Completed: Read {t2m_files[n]}")

    # Open sea ice conc. dataset:   
    ds_sic = xr.open_dataset(sic_files[n])
    if ("valid_time" in ds_sic.coords):
        ds_sic = ds_sic.rename({"valid_time": "time"})
    logging.info(f"Completed: Read {sic_files[n]}") 

    # Open ERA-5 monthly climatological adjustment file:
    t2m_adj_dir = "/dssgfs01/scratch/atb299/ERA5_adjustment/2m_temperature"
    if (n == 1) & (len(ds_t2m.time) == 696): 
       logging.info(f"In Progress: Applying {months[n]} {yr} Leap Year Adjustment.")
       ds_bias = xr.open_dataset(f"{t2m_adj_dir}/2m_temperature-{f'{n+1:02d}'}_LY.nc")
    else:
        logging.info(f"In Progress: Applying {months[n]} {yr} Regular Adjustment.")
        ds_bias = xr.open_dataset(f"{t2m_adj_dir}/2m_temperature-{f'{n+1:02d}'}.nc")

    # Replace time values with input t2m data:
    ds_bias['time'] = ds_t2m['time']

    # Apply climatological correction where sea ice exists (i.e., sic > 0):
    ds_t2m_adj = xr.where(ds_sic['siconc'] > 0, ds_t2m['t2m'] - ds_bias['t2m'], ds_t2m['t2m']).to_dataset(name='t2m')

    # Add attributes of original dataset to adjusted dataset:
    ds_t2m_adj.attrs = ds_t2m.attrs

    logging.info("Completed: Applied monthly climatological correction to ERA-5 2m temperature.")

    # Adjust dtypes for latest ERA-5 data:
    if (yr >= 2024):
        ds_t2m_adj['t2m'] = ds_t2m_adj['t2m'].astype('float')

    # -- Write Climatologically Corrected Data -- #
    # Define original encoding for output netCDF file:
    encoding = get_encoding(ds_t2m, yr)
    # Defining output file path:
    outfile = f"{outdir}/2m_temperature_{yr}-{f'{n+1:02d}'}.nc"
    # Write to netCDF file:
    logging.info(f"In Progress: Producing {months[n]} {yr} climatologically corrected ERA-5 2m temperature.")
    ds_t2m_adj.to_netcdf(outfile, encoding=encoding, unlimited_dims='time')
    logging.info(f"Completed: Produced {months[n]} {yr} climatologically corrected ERA-5 2m temperature in {outfile}")
