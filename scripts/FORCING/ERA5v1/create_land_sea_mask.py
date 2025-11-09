"""
create_land_sea_mask.py

Description: This script creates a binary land-sea mask for ERA-5
atmospheric forcing data.

Created By: Jenny Mecking (jmecki@noc.ac.uk)
Last Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import Dependencies -- #
import xarray as xr

# -- Define Input Arguments -- #
# Path to the input land-sea mask file:
# NOTE: All ERA-5 lsm files checked are the same.
infile = '/dssgfs01/scratch/npd/forcing/ERA5/original/2000/land_sea_mask/land_sea_mask_2000-01.nc'
# Path to the output binary land-sea mask file:
outfile = '/dssgfs01/scratch/npd/forcing/ERA5/original/land_sea_mask_binary_0.nc'

# Minimum land fraction required to label grid box land (note - 0 is all ocean and 1 is all land)
land_min = 0

# -- Create Land-Sea Mask -- # 
ds = xr.open_dataset(infile).squeeze()
# Define land-sea binary mask [0 = land, 1 = ocean]:
ds_out = xr.where(ds['lsm'] > land_min, 0, 1).to_dataset(name='lsm')

# -- Write mask data to file -- #
ds_out.to_netcdf(outfile)
