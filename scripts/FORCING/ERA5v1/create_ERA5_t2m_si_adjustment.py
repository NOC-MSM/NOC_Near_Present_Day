"""
create_ERA5_t2m_si_adjustment.py

Description: This script creates daily ERA-5 climatological adjustment files
for 2m air temperature.

Created By: Adam Blaker (atb299@noc.ac.uk)
Last Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import Dependencies -- #
import numpy as np
import xesmf as xe
import xarray as xr
from datetime import timedelta

# -- Define regrid() -- #
def regrid(ds_in: xr.Dataset, ds_out: xr.Dataset, var_in: xr.DataArray, method:str) -> xr.DataArray:
    """
    Regrid xarray DataArray from one grid to another.
    
    Requires coordinates 'lat' and 'lon' to be present in ds_in and ds_out

    var_in must be on the ds_in grid, and will be interpolated to the grid of ds_out
    using the method specified

    Parameters
    ----------
    ds_in : xarray Dataset
        Dataset defining the input grid in terms of ``lat`` and ``lon``.
    ds_out : xarray Dataset
        Dataset defining the output grid in terms of ``lat`` and ``lon``.
    var_in : xarray DataArray
        DataArray to be regridded from input to target output grid.
    method : str
       regridding method, typically 'bilinear'
       see https://xesmf.readthedocs.io/en/latest/notebooks/Compare_algorithms.html
    
    Returns
    -------
    xarray DataArray
        Regridded DataArray from input to target output grid.
    """
    # -- Verify Inputs -- #
    if not isinstance(ds_in, xr.Dataset):
        raise ValueError("Input dataset must be an xarray Dataset.")
    if not isinstance(ds_out, xr.Dataset):
        raise ValueError("Output dataset must be an xarray Dataset.")
    if not isinstance(var_in, xr.DataArray):
        raise ValueError("Input variable must be an xarray DataArray.")
    if not isinstance(method, str):
        raise ValueError("Method must be a string.")
    
    # -- Regrid DataArray -- #
    # Create global regridder object:
    regridder = xe.Regridder(ds_in, ds_out, method, periodic=True)
    # Regrid DataArray:
    var_out = regridder(var_in)

    return var_out

# -- Reads Monthly t2m Climatology produced by make_clim.py -- #
ds_era5_clim = xr.open_dataset('/dssgfs01/scratch/atb299/ERA5_adjustment/ERA5_SAT_1960-2019_monthly_clim.nc')
ds_jra55_clim = xr.open_dataset('/dssgfs01/scratch/atb299/ERA5_adjustment/JRA55_SAT_1960-2019_monthly_clim.nc')

# -- Regrid JRA55-do t2m to ERA5 grid -- #
t2m_jra55_regrd = regrid(ds_jra55_clim, ds_era5_clim, ds_jra55_clim['ts'], 'bilinear')
t2m_bias = (ds_era5_clim['t2m'] - t2m_jra55_regrd).compute()

# Define dates - uses 13 months by repeating January climatology.
t2m_bias_ext = xr.concat([t2m_bias, t2m_bias[0]], dim='month')
dates = xr.cftime_range(start="0000", periods=13, freq="MS", calendar="noleap")
ds_t2m_bias = t2m_bias_ext.to_dataset(name="t2m").rename({"month":"time"})
ds_t2m_bias["time"] = dates

# -- Write monthly t2m climatological bias to netCDF -- #
ds_t2m_bias.isel(time=slice(0,12)).to_netcdf('/dssgfs01/scratch/atb299/ERA5_adjustment/SAT_clim_diff_E-J.nc')

# Resample to daily and roll the calendar and data so values are centre of the month
ds_t2m_bias_daily = ds_t2m_bias.resample(time="1D").interpolate("linear")
ds_t2m_bias_daily["time"] = ds_t2m_bias_daily.get_index("time") + timedelta(days=15)
ds_t2m_bias_daily['t2m'] = ds_t2m_bias_daily['t2m'].roll(time=15)

# -- Write daily t2m climatological bias to netCDF -- #
ds_t2m_bias_daily.to_netcdf('/dssgfs01/scratch/atb299/ERA5_adjustment/SAT_clim_diff_daily_E-J.nc')

# Resample to create hourly data in monthly files
mond = [0,31,28,31,30,31,30,31,31,30,31,30,31]

# Iterate over months:
for i in range(12):
    # Non-Leap Year Day Indexes:
    tmin=np.sum(mond[:i+1]).astype(int)
    tmax=(np.sum(mond[:i+2])+1).astype(int)
    print("Days of year:", tmin, tmax)
    # Calculate hourly t2m data for each month:
    t2m_bias_hourly = ds_t2m_bias_daily.isel(time=slice(tmin,tmax)).resample(time="1h").interpolate("linear").isel(time=slice(None,-1))
    # Write hourly t2m data to netCDF:
    t2m_bias_hourly.to_netcdf('/dssgfs01/scratch/atb299/ERA5_adjustment/2m_temperature/2m_temperature-'+f"{i+1:02d}"+'.nc')
    print(f"Month {i} shape is: {t2m_bias_hourly['t2m'].shape}")
    
    if i == 1:
        # Calculate hourly t2m data for February in leap years:
        t2m_bias_hourly = ds_t2m_bias_daily.isel(time=slice(tmin,tmax+1)).resample(time="1h").interpolate("linear").isel(time=slice(None,-1))
        # Write hourly t2m data to netCDF:
        t2m_bias_hourly.to_netcdf('/dssgfs01/scratch/atb299/ERA5_adjustment/2m_temperature/2m_temperature-'+f"{i+1:02d}"+'_LY.nc')
        print(f"Month {i} leap year shape is: {t2m_bias_hourly['t2m'].shape}")
