"""
utils.py

Utility functions for generating initial conditions for NEMO ocean 
general circulation models using World Ocean Atlas climate norms.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
Created On: 25/11/2024
"""
# -- Import Python Packages -- #
import xoak
from tqdm import tqdm
import cf_xarray
import numpy as np
import xarray as xr 

# -- Determine indexes of nearest non-NaN grid points to NaN grid points -- #
def get_ffill_indexes(ds:xr.Dataset, var:str = 't_an') -> xr.Dataset:
    """
    Find indexes (y,x) of nearest non-NaN grid points to NaN grid points
    to flood fill NaN values at each depth level.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing 1D longitude (lon) and latitude (lat) arrays
        and variable (e.g. temperature / salinity) containing NaNs.
    var : str
        Name of variable to be filled.
    """
    # -- Verify input arguments -- #
    if isinstance(ds, xr.Dataset) is False:
        raise TypeError("ds must be an xarray Dataset.")
    if 'lon' not in ds.coords or 'lat' not in ds.coords:
        raise ValueError("ds must contain longitude (lon) and latitude (lat) coordinates.")
    if 'depth' not in ds.dims:
        raise ValueError("ds must contain depth dimension.")
    if isinstance(var, str) is False:
        raise TypeError("var must be a string.")
    if var not in ds:
        raise ValueError(f"{var} not found in ds.")

    # -- Define 2D (lat, lon) grid -- #
    # Define 2-dimensional longitude and latitude arrays:
    lon, lat = np.meshgrid(ds.lon.values, ds.lat.values)

    # Create a Dataset with 2D lon and lat arrays:
    ds_3D_grid = xr.Dataset(coords={'lat': (('y', 'x'), lat), 'lon': (('y', 'x'), lon)})
    # Add x/y_index coordinates:
    for key, value in ds_3D_grid.sizes.items():
        ds_3D_grid[f"{key}_index"] = xr.DataArray(range(value), dims=key)
    ds_3D_grid = ds_3D_grid.cf.guess_coord_axis()

    # Expand dimensions to (depth,y,x):
    ds_3D_grid['x_index'] = ds_3D_grid['x_index'].expand_dims(dim={"depth":ds.sizes["depth"], "y": ds_3D_grid.sizes['y']}).transpose('depth', 'y', 'x')
    ds_3D_grid['y_index'] = ds_3D_grid['y_index'].expand_dims(dim={"depth":ds.sizes["depth"], "x": ds_3D_grid.sizes['x']}).transpose('depth', 'y', 'x')

    # Define 3D index arrays in terms of (depth,y,x) using deep copy:
    x_3D_index = xr.DataArray(ds_3D_grid['x_index'].values.copy(), dims=('depth', 'y', 'x'))
    y_3D_index = xr.DataArray(ds_3D_grid['y_index'].values.copy(), dims=('depth', 'y', 'x'))

    # Iterate over WOA depth levels:
    for nd in tqdm(range(ds.sizes['depth'])):
        # Subset the WOA temperature array to 2D (y,x):
        t_an_2D = ds[var].squeeze().isel(depth=nd).values
        # Create a mask for land and sea points on WOA grid:
        mask_land = np.where(np.isnan(t_an_2D))
        mask_sea = np.where(~np.isnan(t_an_2D))
        # Extract 1D longitude and latitude arrays for sea points:
        lon_sea = lon[mask_sea].flatten()
        lat_sea = lat[mask_sea].flatten()

        # Create a Dataset with the longitude and latitude arrays for sea points:
        ds_1D_grid = xr.Dataset(coords={'lat': (('xs'), lat_sea), 'lon': (('xs'), lon_sea)})
        # Add 1D index coordinate xs:
        for key, value in ds_1D_grid.sizes.items():
            ds_1D_grid[f"{key}_index"] = xr.DataArray(range(value), dims=key)
        ds_1D_grid = ds_1D_grid.cf.guess_coord_axis()

        # Set index using xoak using the py2sindex based on s2 geometry:
        ds_1D_grid.xoak.set_index(('lat', 'lon'), 's2point')

        # Extract 1D longitude and latitude arrays for land points:
        lon_land = lon[mask_land].flatten()
        lat_land = lat[mask_land].flatten()

        # Get the 1D index coordinates of the nearest sea point to all land points:
        ds_index = ds_1D_grid.xoak.sel(lat=xr.DataArray(lat_land), lon=xr.DataArray(lon_land))

        # Transform 1D index coordinates in (xs) to 2D index coordinates in (y,x):
        x_index = ds_3D_grid.x_index.isel(depth=nd).values[mask_sea][ds_index.xs_index.values]
        y_index = ds_3D_grid.y_index.isel(depth=nd).values[mask_sea][ds_index.xs_index.values]

        # Update the 2D index coordinates (y,x) of land points with nearest sea point (y,x):
        x_3D_index.isel(depth=nd).values[mask_land] = ds_3D_grid.x_index.isel(depth=nd).values[y_index, x_index]
        y_3D_index.isel(depth=nd).values[mask_land] = ds_3D_grid.y_index.isel(depth=nd).values[y_index, x_index]

        ds_3D_index = xr.Dataset({'x_index': x_3D_index, 'y_index': y_3D_index})

    return ds_3D_index