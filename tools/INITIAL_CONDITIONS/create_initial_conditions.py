"""
create_initial_conditions.py

Workflow for generating initial conditions for NEMO ocean general
circulation models by flood filling and interpolating the temperature
and salinity fields provided by World Ocean Atlas climate norms.

== Steps are as follows ==

1. Load World Ocean Atlas (WOA23) temperature and salinity data.
2. Compute or read flood fill indexes for NaN WOA23 land grid cells.
3. Merge full-depth and upper 1500m WOA23 temperature and salinity fields.
4. Flood fill NaN values in merged temperature and salinity fields.
5. Vertically interpolate WOA23 temperature and salinity fields onto NEMO eORCA grid.
6. Calculate TEOS-10 Conservative Temperature and Absolute Salinity.
7. Regrid TEOS-10 Conservative Temperature and Absolute Salinity onto NEMO eORCA grid.
8. Save initial conditions to NetCDF file.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
Created On: 25/11/2024
"""
# ============================================ DEFINE USER INPUTS ============================================ #
# Define filepaths for WOA23 and NEMO model data:
working_dir = "/dssgfs01/working/otooth/NOC_NPD/20241017/NOC_Near_Present_Day/tools/INITIAL_CONDITIONS/data"
woa23_wa_temp_filepath = f"{working_dir}/woa23_decav_t13_04.nc"
woa23_wa_sal_filepath = f"{working_dir}/woa23_decav_s13_04.nc"
nemo_grid_filepath = "/dssgfs01/working/atb299/NEMO_cfgs/eORCA025_RK3/mesh_mask.nc"

# Define climate norms decade for WOA23 temperature and salinity data:
woa23_decav = "decav71A0"

# -- Define Boolean Flags -- #
# Set to True to calculate flood fill indexes:
calc_ffill_indexes = False
# Set to True to save flood fill indexes to netCDF file:
save_ffill_indexes = False
# Define filepath to read flood fill indexes from netCDF file - otherwise set to None:
ffill_indexes_filepath = f"{working_dir}/woa23_ffill_indexes.nc"

# -- Define output netCDF file parameters -- #
# Define output chunk sizes:
out_chunks = (1, 5, 64, 64)
# Define output file path:
out_filepath = f"{working_dir}/woa23_decav71A0_TS_TEOS10_eORCA025.nc"

# ============================================ END OF USER INPUTS ============================================ #

# -- Import Python Packages -- #
import sys
import gsw
import logging
import xesmf as xe
import numpy as np
import xarray as xr

# -- Import User Defined Functions -- #
from utils import get_ffill_indexes

# -- Step 0: Set up Logging -- #
logger = logging.getLogger(__name__)
logging.basicConfig(
        stream=sys.stdout,
        format="~ Initial Conditions | %(levelname)10s | %(asctime)s | %(message)s",
        level=logging.INFO,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

# -- Step 1: Load World Ocean Atlas Data -- #
logging.info("Step 1: Reading full-depth World Ocean Atlas Data.")
# Import WOA23 winter-average gridded temperature and salinity data:
ds_full_temp = xr.open_dataset(woa23_wa_temp_filepath, decode_times=False)['t_an'].to_dataset()
ds_full_sal = xr.open_dataset(woa23_wa_sal_filepath, decode_times=False)['s_an'].to_dataset()

# Define 2-dimensional longitude and latitude arrays:
lon, lat = np.meshgrid(ds_full_temp.lon.values, ds_full_temp.lat.values)
# Construct Dataset with 2D longitude and latitude arrays:
ds_full_TS = xr.Dataset(
    data_vars={'t_an': (('depth', 'y','x'), ds_full_temp.t_an.squeeze().values),
               's_an':(('depth', 'y','x'), ds_full_sal.s_an.squeeze().values)},
    coords={'depth':(('depth'), ds_full_temp.depth.data),
            'lat': (('y','x'), lat),
            'lon': (('y','x'), lon)}
            )

# Expand time dimension to 12 months to merge WOA23 temperature and salinity datasets:
ds_full_TS = ds_full_TS.expand_dims(dim={"time":12})

# -- Step 2 (Optional): Get flood-fill indexes for NaN land grid cells -- #
if calc_ffill_indexes:
    logging.info("Step 2: Calculating 3D indexes to flood fill World Ocean Atlas Data on land.")
    # Get the 3D indexes to flood fill temperature and salinity NaN values on land:
    ds_3D_index = get_ffill_indexes(ds=ds_full_temp, var='t_an')

    if save_ffill_indexes:
        # Write the 3D index arrays Dataset to netCDF file:
        out_filepath = f"{working_dir}/woa23_ffill_indexes.nc"
        ds_3D_index.to_netcdf(out_filepath)

else:
    logging.info("Step 2: Reading 3D indexes to flood fill World Ocean Atlas Data on land from netCDF file.")
    if ffill_indexes_filepath is not None:
        # Read the 3D index arrays Dataset from netCDF file:
        ds_3D_index = xr.open_dataset(ffill_indexes_filepath)

# -- Step 3: Merge WOA23 temperature and salinity fields -- #
logging.info("Step 3: Merge full-depth and upper 1500m World Ocean Atlas data.")
# Import Upper 1500m WOA gridded temperature data:
# Define file paths for temperature and salinity data:
woa23_1500m_temp_filepath=f"{working_dir}/woa23_{woa23_decav}_t*_04.nc"
woa23_1500m_sal_filepath=f"{working_dir}/woa23_{woa23_decav}_s*_04.nc"

# Open WOA23 monthly temperature and salinity climatology files as multi-file dataset:
ds_u1500_temp = xr.open_mfdataset(woa23_1500m_temp_filepath, decode_times=False, preprocess=lambda ds: ds['t_an'])
ds_u1500_sal = xr.open_mfdataset(woa23_1500m_sal_filepath, decode_times=False, preprocess=lambda ds: ds['s_an'])

# Define 2-dimensional longitude and latitude arrays (y,x):
lon, lat = np.meshgrid(ds_u1500_temp.lon.values, ds_u1500_temp.lat.values)
# Construct Dataset with 2D longitude and latitude arrays & load into memory:
ds_u1500_grid = xr.Dataset(
    data_vars={'t_an': (('time', 'depth', 'y', 'x'), ds_u1500_temp.t_an.data),
               's_an': (('time', 'depth', 'y', 'x'), ds_u1500_sal.s_an.data)},
    coords={'time':(('time'), ds_u1500_temp.time.data),
            'depth':(('depth'), ds_u1500_temp.depth.data),
            'lat': (('y', 'x'), lat),
            'lon': (('y', 'x'), lon)}
            ).load()

# Store time values:
time = ds_u1500_grid.time.values
# Remove used datasets from memory:
del ds_u1500_temp, ds_u1500_sal, ds_full_temp, ds_full_sal

# Replace full-depth temperature and salinity data with upper 1500m (upper 57 grid cells) data:
ds_full_TS['t_an'] = ds_full_TS['t_an'].copy()
ds_full_TS['s_an'] = ds_full_TS['s_an'].copy()
ds_full_TS['s_an'][:,:57,:,:] = ds_u1500_grid['s_an']
ds_full_TS['t_an'][:,:57,:,:] = ds_u1500_grid['t_an']

# Remove used datasets from memory:
del ds_u1500_grid

# -- Step 4: Flood Fill Temperature and Salinity Fields -- #
logging.info("Step 4: Flood fill land in merged World Ocean Atlas temperature and salinity fields.")
# Flood fill land NaN values in temperature and salinity data with nearest sea point values:
# 1. Full-depth Temperature:
ds_full_TS['t_an'] = ds_full_TS['t_an'].isel(y=ds_3D_index.y_index, x=ds_3D_index.x_index)
# 2. Full-depth Salinity:
ds_full_TS['s_an'] = ds_full_TS['s_an'].isel(y=ds_3D_index.y_index, x=ds_3D_index.x_index)

# -- Step 5: Vertical Interpolation of Temperature and Salinity Fields -- #
logging.info("Step 5: Vertically interpolate merged, flood filled World Ocean Atlas data onto eORCA grid.")
# Read target NEMO eORCA grid mesh_mask file:
ds_mesh_mask = xr.open_dataset(nemo_grid_filepath)

# Vertical interpolation of temperature and salinity data to NEMO eORCA vertical grid:
# NOTE: Interpolation fill_value set to 'extrapolate' to avoid NaN values.
ds_full_TS_vinterp = ds_full_TS.interp(depth=ds_mesh_mask['gdept_1d'].squeeze().values,
                                       method='linear',
                                       kwargs={'fill_value': 'extrapolate'}
                                        )

# Replace depth coord and convert data types to float32:
for var in ds_full_TS_vinterp.data_vars:
    ds_full_TS_vinterp[var] = ds_full_TS_vinterp[var].astype('float32')
# Rename variable depth to deptht:
ds_full_TS_vinterp = ds_full_TS_vinterp.rename({'depth': 'deptht'})

# Remove used datasets from memory:
del ds_full_TS

# -- Step 6: Calculate TEOS-10 Conservative Temperature and Absolute Salinity -- #
logging.info("Step 6: Convert World Ocean Atlas temperature and salinity to TEOS-10.")
# Compute TEOS-10 Absolute Salinity - note latitudes cannot south of than -86N otherwise NaN, so we fill with -86N:
ds_full_TS_vinterp['s_abs']=gsw.SA_from_SP(SP=ds_full_TS_vinterp['s_an'],
                                           p=ds_full_TS_vinterp['deptht'],
                                           lon=ds_full_TS_vinterp['lon'],
                                           lat=ds_full_TS_vinterp['lat'].where(cond=ds_full_TS_vinterp['lat'] >= -86, other=-86))

# Compute TEOS-10 Conservative Temperature:
ds_full_TS_vinterp['t_con'] = gsw.CT_from_t(SA=ds_full_TS_vinterp['s_abs'], t=ds_full_TS_vinterp['t_an'], p=ds_full_TS_vinterp['deptht'])

# -- Step 7: Regrid to NEMO eORCA grid -- #
logging.info("Step 7: Regrid TEOS-10 World Ocean Atlas temperature and salinity onto eORCA grid.")
# Define xesmf regrid function from WOA23 grid to NEMO eORCA grid:
regridder = xe.Regridder(ds_full_TS_vinterp,
                         ds_mesh_mask.rename({'nav_lon':'lon','nav_lat':'lat'}),
                         'bilinear',
                         periodic=True
                         )

# Regrid WOA23 temperature and salinity fields to NEMO eORCA grid:
ds_TS_out = regridder(ds_full_TS_vinterp['s_abs']).to_dataset(name='salinity')
ds_TS_out['temperature'] = regridder(ds_full_TS_vinterp['t_con'])

# -- Step 8: Save Initial Conditions to NetCDF File -- #
logging.info("Step 8: Writing regrided TEOS-10 World Ocean Atlas temperature and salinity to netCDF file.")
# Add time coordinate with values from WOA23 temperature and salinity data:
ds_TS_out = ds_TS_out.assign_coords({'time':(('time'), time)})
# Add NEMO eORCA grid coordinates:
ds_TS_out['nav_lat'] = ds_mesh_mask['nav_lat']
ds_TS_out['nav_lon'] = ds_mesh_mask['nav_lon']
# Convert data types to float32:
for var in ds_TS_out.data_vars:
    ds_TS_out[var] = ds_TS_out[var].astype('float32')

ds_TS_out.to_netcdf(out_filepath, encoding={"temperature": {"chunksizes": out_chunks}, "salinity": {"chunksizes": out_chunks}})
logging.info("✔ Created WOA23 initial conditions for eORCA grid successfully ✔")

# ============================================ END OF SCRIPT ============================================ #