"""
create_sss_restoring_climatology.py

Workflow for generating sea surface salinity restoring climatology for NEMO ocean
general circulation models by flood filling and interpolating the upper 10m mean
salinity field provided by World Ocean Atlas climate norms.

== Steps are as follows ==

1. Load World Ocean Atlas (WOA23) upper 1500m salinity data.
2. Compute or read flood fill indexes for NaN WOA23 land grid cells.
3. Flood fill NaN values in merged salinity fields.
4. Compute depth-average of upper 10m salinity field.
5. Calculate TEOS-10 upper 10m average Absolute Salinity.
7. Regrid TEOS-10 upper 10m average Absolute Salinity onto NEMO eORCA grid.
8. Save SSS restoring climatology to NetCDF file.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
Created On: 06/12/2024
"""
# ============================================ DEFINE USER INPUTS ============================================ #
# Define file directory for WOA23 data:
working_dir = "/dssgfs01/working/otooth/NOC_NPD/20241017/NOC_Near_Present_Day/tools/INITIAL_CONDITIONS/data"

# Define NEMO grid file path:
nemo_grid_filepath = "/dssgfs01/working/atb299/NEMO_cfgs/NPD_eORCA1_v4.2/domain_cfg.nc"
# nemo_grid_filepath = "/dssgfs01/working/atb299/NEMO_cfgs/eORCA025_RK3/mesh_mask.nc"
# nemo_grid_filepath = "/dssgfs01/working/atb299/NEMO_cfgs/NPD_eORCA12_v4.2/domain_cfg.nc"

# Define climate norms decade for WOA23 salinity data:
woa23_decav = "decav91C0"
# Define upper 1500m filepaths:
woa23_1500m_sal_filepath=f"{working_dir}/woa23_{woa23_decav}_s*_04.nc"

# -- Define Boolean Flags -- #
# Set to True to calculate flood fill indexes:
calc_ffill_indexes = False
# Set to True to save flood fill indexes to netCDF file:
save_ffill_indexes = False
# Define filepath to read flood fill indexes from netCDF file - otherwise set to None:
ffill_indexes_filepath = f"{working_dir}/woa23_ffill_indexes.nc"

# -- Define output netCDF file parameters -- #
# Define output chunk sizes:
out_chunks = (1, 32, 32)
# Define output file path:
out_filepath = f"{working_dir}/woa23_{woa23_decav}_sss_TEOS10_eORCA1.nc"
# out_filepath = f"{working_dir}/woa23_{woa23_decav}_sss_TEOS10_eORCA025.nc"
# out_filepath = f"{working_dir}/woa23_{woa23_decav}_sss_TEOS10_eORCA12.nc"

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
def setup_logging():
    """
    Set up logging configuration to output to stdout.
    """
    logger = logging.getLogger(__name__)
    logging.basicConfig(
            stream=sys.stdout,
            format="~ SSS Restoring | %(levelname)10s | %(asctime)s | %(message)s",
            level=logging.INFO,
            datefmt="%Y-%m-%d %H:%M:%S",
        )

# -- Step 1: Load WOA23 upper 1500m salinity fields -- #
def load_woa23_data(woa23_1500m_sal_filepath):
    """
    Load World Ocean Atlas 2023 upper 1500m salinity data.
    """
    logging.info("Step 1: Reading upper 1500m World Ocean Atlas salinity data.")
    # Open WOA23 monthly upper 10m [0,5,10] salinity climatology files as multi-file dataset:
    ds_u10_sal = xr.open_mfdataset(woa23_1500m_sal_filepath, decode_times=False, preprocess=lambda ds: ds['s_an'].isel(depth=slice(0, 3)))

    # Define 2-dimensional longitude and latitude arrays (y,x):
    lon, lat = np.meshgrid(ds_u10_sal.lon.values, ds_u10_sal.lat.values)
    # Construct Dataset with 2D longitude and latitude arrays & load into memory:
    ds_u10_grid = xr.Dataset(
        data_vars={'s_an': (('time', 'depth', 'y', 'x'), ds_u10_sal.s_an.data)},
        coords={'time':(('time'), ds_u10_sal.time.data),
                'depth':(('depth'), ds_u10_sal.depth.data),
                'lat': (('y', 'x'), lat),
                'lon': (('y', 'x'), lon)}
                ).load()

    # Store time values:
    time = ds_u10_grid.time.values

    return ds_u10_grid, time

# -- Step 2 (Optional): Get flood-fill indexes for NaN land grid cells -- #
def prepare_ffill_indexes(calc_ffill_indexes, ds_u10_sal, save_ffill_indexes, ffill_indexes_filepath):
    """
    Calculate or read 3D indexes to flood fill World Ocean Atlas 2023 data on land.
    """
    if calc_ffill_indexes:
        logging.info("Step 2: Calculating 3D indexes to flood fill World Ocean Atlas Data on land.")
        # Get the 3D indexes to flood fill temperature and salinity NaN values on land:
        ds_3D_index = get_ffill_indexes(ds=ds_u10_sal.isel(time=0), var='s_an')

        if save_ffill_indexes:
            # Write the 3D index arrays Dataset to netCDF file:
            ds_3D_index.to_netcdf(ffill_indexes_filepath)

    else:
        logging.info("Step 2: Reading 3D indexes to flood fill World Ocean Atlas Data on land from netCDF file.")
        if ffill_indexes_filepath is not None:
            # Read the 3D index arrays Dataset from netCDF file:
            ds_3D_index = xr.open_dataset(ffill_indexes_filepath)

    return ds_3D_index

# -- Step 3: Flood Fill Salinity Fields -- #
def ffill_woa23_data(ds_u10_grid, ds_3D_index):
    """
    Flood fill NaN values in merged World Ocean Atlas 2023 temperature and salinity fields.
    """
    logging.info("Step 3: Flood fill land in merged World Ocean Atlas temperature and salinity fields.")
    # Flood fill land NaN values in salinity data with nearest sea point values:
    # Upper 1500m Salinity:
    ds_u10_grid['s_an'] = ds_u10_grid['s_an'].isel(y=ds_3D_index.y_index.isel(depth=slice(0, 3)), x=ds_3D_index.x_index.isel(depth=slice(0, 3)))

    return ds_u10_grid

# -- Step 4: Compute Depth-Average Salinity Field -- #
def calc_depth_avg_salinity(ds_u10_grid):
    """
    Calculate upper 10m depth-average of flood filled World Ocean Atlas 2023 salinity.
    """
    logging.info("Step 4: Calculate upper 10m depth-average of flood filled World Ocean Atlas salinity.")
    # Read target NEMO eORCA grid mesh_mask file:
    ds_u10_grid['s_abs'] = ds_u10_grid['s_an'].mean(dim='depth')

    # Replace depth coord and convert data types to float32:
    for var in ds_u10_grid.data_vars:
        ds_u10_grid[var] = ds_u10_grid[var].astype('float32')
    # Rename variable depth to deptht:
    ds_u10_grid = ds_u10_grid.rename({'depth': 'deptht'})

    return ds_u10_grid

# -- Step 5: Calculate TEOS-10 Conservative Temperature and Absolute Salinity -- #
def convert_to_TEOS10(ds_u10_grid):
    """
    Convert World Ocean Atlas 2023 upper 10m average salinity to TEOS-10.
    """
    logging.info("Step 5: Convert World Ocean Atlas upper 10m average salinity to TEOS-10.")
    # Compute TEOS-10 Absolute Salinity - note latitudes cannot south of than -86N otherwise NaN, so we fill with -86N:
    ds_u10_grid['s_abs']=gsw.SA_from_SP(SP=ds_u10_grid['s_abs'],
                                        p=5,
                                        lon=ds_u10_grid['lon'],
                                        lat=ds_u10_grid['lat'].where(cond=ds_u10_grid['lat'] >= -86, other=-86))
    return ds_u10_grid

# -- Step 6: Regrid to NEMO eORCA grid -- #
def regrid_woa23_data(ds_u10_grid, nemo_grid_filepath):
    """
    Regrid TEOS-10 World Ocean Atlas upper 10m average salinity onto NEMO eORCA grid.
    """
    logging.info("Step 6: Regrid TEOS-10 World Ocean Atlas upper 10m average salinity onto eORCA grid.")
    # Read target NEMO eORCA grid mesh_mask file:
    ds_mesh_mask = xr.open_dataset(nemo_grid_filepath)
    ds_mesh_mask['nav_lon'] = ds_mesh_mask['glamt'].squeeze()
    ds_mesh_mask['nav_lat'] = ds_mesh_mask['gphit'].squeeze()

    # Define xesmf regrid function from WOA23 grid to NEMO eORCA grid:
    regridder = xe.Regridder(ds_u10_grid,
                            ds_mesh_mask.rename({'nav_lon':'lon','nav_lat':'lat'}),
                            'bilinear',
                            periodic=True
                            )

    # Regrid WOA23 temperature and salinity fields to NEMO eORCA grid:
    ds_sss_out = regridder(ds_u10_grid['s_abs']).to_dataset(name='salinity')

    return ds_sss_out, ds_mesh_mask

# -- Step 7: Save Initial Conditions to NetCDF File -- #
def write_regridded_woa23_data_to_dataset(ds_sss_out, ds_mesh_mask, time, out_filepath, out_chunks):
    """
    Save regrided TEOS-10 World Ocean Atlas upper 10m average salinity to netCDF file.
    """
    logging.info("Step 7: Writing regrided TEOS-10 World Ocean Atlas upper 10m average salinity to netCDF file.")
    # Add time coordinate with values from WOA23 temperature and salinity data:
    ds_sss_out = ds_sss_out.assign_coords({'time':(('time'), time), 'nav_lat':(('y', 'x'), ds_mesh_mask['nav_lat'].data), 'nav_lon':(('y', 'x'), ds_mesh_mask['nav_lon'].data)})
    # Convert data types to float32:
    for var in ds_sss_out.data_vars:
        ds_sss_out[var] = ds_sss_out[var].astype('float32')

    # Write output to netCDF file:
    ds_sss_out.to_netcdf(out_filepath, encoding={"salinity": {"chunksizes":out_chunks, "zlib":True, "complevel":4}}, unlimited_dims="time")

# ============================================ MAIN WORKFLOW ============================================ #
if __name__ == "__main__":
    # Step 0: Set up logging configuration:
    setup_logging()
    # Step 1: Load World Ocean Atlas 2023 upper 1500m salinity data:
    ds_u10_grid, time = load_woa23_data(woa23_1500m_sal_filepath)
    # Step 2: Get flood-fill indexes for NaN land grid cells:
    ds_3D_index = prepare_ffill_indexes(calc_ffill_indexes, ds_u10_grid['s_an'], save_ffill_indexes, ffill_indexes_filepath)
    # Step 3: Flood fill NaN values in merged World Ocean Atlas 2023 salinity fields:
    ds_u10_grid = ffill_woa23_data(ds_u10_grid, ds_3D_index)
    # Step 4: Compute depth-average of upper 10m salinity field:
    ds_u10_grid = calc_depth_avg_salinity(ds_u10_grid)
    # Step 5: Calculate TEOS-10 upper 10m average Absolute Salinity:
    ds_u10_grid = convert_to_TEOS10(ds_u10_grid)
    # Step 6: Regrid TEOS-10 upper 10m average Absolute Salinity onto NEMO eORCA grid:
    ds_sss_out, ds_mesh_mask = regrid_woa23_data(ds_u10_grid, nemo_grid_filepath)
    # Step 7: Save SSS restoring climatology to NetCDF file:
    write_regridded_woa23_data_to_dataset(ds_sss_out, ds_mesh_mask, time, out_filepath, out_chunks)

    logging.info("✔ Created WOA23 SSS restoring climatology for eORCA grid successfully ✔")

# ============================================ MAIN WORKFLOW ============================================ #