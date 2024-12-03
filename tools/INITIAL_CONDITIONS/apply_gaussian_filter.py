"""
apply_gaussian_filter.py

Workflow for applying a Gaussian smoothing filter to WOA23 initial conditions
regridded for NEMO ocean general circulation models. Gaussian filtering is used
to smooth temperature and salinity gradients arising from nearest-neighbour
interpolation in create_initial_conditions.py.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
Created On: 03/12/2024
"""
# ============================================ DEFINE USER INPUTS ============================================ #
# Define filepath to WOA23 flood-filled interpolated initial conditions:
filepath = "/dssgfs01/working/otooth/NEMO_cfgs/NPD_eORCA025_v4.2/woa23_decav71A0_TS_TEOS10_eORCA025.nc"

# Define x-y subdomain limits to apply Gaussian smoothing:
xmin, xmax = 1242, 1270
ymin, ymax = 860, 880

# Define sigma value for Gaussian filter:
n_pixel_sigma = 3

# -- Define output netCDF file parameters -- #
# Define output chunk sizes:
out_chunks = (1, 5, 32, 32)

# ============================================ END OF USER INPUTS ============================================ #

# -- Import Python Packages -- #
import sys
import logging
import scipy
import xarray as xr
from tqdm import tqdm

# -- Step 0: Set up Logging -- #
logger = logging.getLogger(__name__)
logging.basicConfig(
        stream=sys.stdout,
        format="~ Gaussian Filter | %(levelname)10s | %(asctime)s | %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

# -- Step 1: Load World Ocean Atlas Data -- #
logging.info("Step 1: Reading regridded World Ocean Atlas data.")
# Import regridded WOA23 initial conditions:
ds = xr.open_dataset(filepath)

# -- Step 2: Preparing temperature and salinity data -- #
logging.info("Step 2: Preparing regridded temperature and salinity data.")
# Define number of depth levels:
n_depth = ds.deptht.size
# Define number of time steps:
n_time = ds.time.size

# Extract regridded WOA23 temperature and salinity fields:
salinity = ds.salinity.values.copy()
temperature = ds.temperature.values.copy()

# -- Step 3: Apply Gaussian filter to temperature and salinity data -- #
logging.info("Step 3: Apply Gaussian filter to subdomain of temperature and salinity data.")
# Iterate over each time step:
for nt in tqdm(range(n_time)):
# Apply a Gaussian filter to the T-S fields at each depth level:
    for nd in range(n_depth):
        # Replace the values of the absolute salinity field with the filtered values:
        salinity[nt, nd, ymin:ymax, xmin:xmax] = scipy.ndimage.gaussian_filter(salinity[nt, nd, ymin:ymax, xmin:xmax], sigma=n_pixel_sigma, order=0)

        # Replace the values of the conservative temperature field with the filtered values:
        temperature[nt, nd, ymin:ymax, xmin:xmax] = scipy.ndimage.gaussian_filter(temperature[nt, nd, ymin:ymax, xmin:xmax], sigma=n_pixel_sigma, order=0)

# -- Step 4: Apply Gaussian filter to temperature and salinity data -- #
logging.info("Step 4: Writing Gaussian filtered temperature and salinity data to netCDF file.")
# Define the Gaussian filtered absolute salinity and conservative temperature WOA23 fields as DataArrays:
salinity = xr.DataArray(salinity, coords=[ds.time, ds.deptht, ds.y, ds.x], dims=['time', 'deptht', 'y', 'x'])
temperature = xr.DataArray(temperature, coords=[ds.time, ds.deptht, ds.y, ds.x], dims=['time', 'deptht', 'y', 'x'])

# Write the Gaussian filtered absolute salinity and conservative temperature WOA23 fields to a new netCDF file:
ds_out = xr.Dataset({'salinity': salinity, 'temperature': temperature})
ds_out.to_netcdf(filepath.replace('.nc', '_gauss_smoothed.nc'), encoding={"temperature": {"chunksizes": out_chunks}, "salinity": {"chunksizes": out_chunks}}, unlimited_dims="time")

logging.info("✔ Applied Gaussian smoothing filter to regridded WOA23 data successfully ✔")
