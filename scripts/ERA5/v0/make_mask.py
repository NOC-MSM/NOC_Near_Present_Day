import xarray as xr
import numpy as np

# File on which to base the mask (note - all ERA5 lsm files checked the are the same):
infile   = '/dssgfs01/scratch/npd/forcing/ERA5/original/2000/land_sea_mask/land_sea_mask_2000-01.nc'
outfile  = '/dssgfs01/scratch/npd/forcing/ERA5/original/land_sea_mask_binary_0_test.nc'
# Minimum land fraction required to label grid box land (note - 0 is all ocean and 1 is all land)
land_min = 0

# Read in data: 
ds = xr.open_dataset(infile)

# Collapse to only have 1 time step:
ds_new = ds.mean(dim='time')

# Change to a binary mask:
ds_new.lsm.values[np.where(ds_new.lsm.values>land_min)] = 1
ds_new.lsm.values[np.where(ds_new.lsm.values!=1)] = 0

# Make land 0 and ocean 1:
ds_new.lsm.values = 1- ds_new.lsm.values

# Save data to file:
ds_new.to_netcdf(outfile)
