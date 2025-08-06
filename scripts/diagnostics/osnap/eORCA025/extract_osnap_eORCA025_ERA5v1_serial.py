# ---------------------------------------------------
# extract_osnap_eORCA025_ERA5_serial.py
# 
# Description:
# Script to extract the OSNAP trans-basin section from
# the eORCA025 ERA-5 version 1 Near-Present Day
# simulation.
# 
# Created By:
# Ollie Tooth (oliver.tooth@noc.ac.uk) 
#
# ----------------------------------------------------
# === Import Python packages === #
import gsw
import glob
import numpy as np
import xarray as xr
from nemo_cookbook import extract_section, compute_section_moc_tracer

# === Define function to extract OSNAP section and compute MOC === #
def extract_osnap_moc(domain_path: str,
                      T_paths: dict[str, list[str]],
                      U_paths: dict[str, list[str]],
                      V_paths: dict[str, list[str]],
                      out_file: str,
                      ) -> None:
    """
    Extract the Overturning in the Subpolar North Atlantic Program (OSNAP)
    section from the eORCA025 ERA5 v1 simulation and compute the meridional
    overturning stream function in potential density coordinates.

    Parameters
    ----------
    domain_path : str
        Path to the NEMO domain configuration file.
    T_paths : str
        Regular expression defining file paths to the NEMO grid T files.
    U_paths : str
        Regular expression defining to the NEMO grid U files.
    V_paths : str
        Regular expression defining to the NEMO grid V files.
    out_file : str
        Path to the output NetCDF file where the OSNAP section will be saved.
    """
    # === Initialise Dask Local Cluster === #
    print("=== Extracting eORCA025 ERA-5 v1 OSNAP MOC without Dask ===")

    # === Prepare OSNAP Coordinates === #
    # Define S3 url to OSNAP gridded observational data in JASMIN Object Store:
    url = "https://noc-msm-o.s3-ext.jc.rl.ac.uk/ocean-obs/OSNAP/OSNAP_gridded_2014_2020/"
    ds_osnap = xr.open_zarr(url, zarr_format=3, consolidated=True)

    # Define observation coordinates defining the OSNAP array:
    # Note: An final coordinate is added to ensure the Scottish shelf is included in the section.
    osnap_lons = np.concatenate([ds_osnap['LONGITUDE'].values[::3], np.array([-4.5])])
    osnap_lats = np.concatenate([ds_osnap['LATITUDE'].values[::3], np.array([56.5])])

    # === Extract OSNAP Section === #
    # Extract the OSNAP section from the NEMO model data:
    ds_osnap = extract_section(section_lon=osnap_lons,
                                section_lat=osnap_lats,
                                domain_path=domain_path,
                                T_paths=T_paths,
                                U_paths=U_paths,
                                V_paths=V_paths,
                                var_map={'temp':'thetao_con', 'sal':'so_abs'},
                                uv_eiv=True,
                                log=True,
                                )

    # === Compute OSNAP AMOC in density-coords === #
    # Calculate potential density referenced to the sea surface using the Gibbs Sea Water Toolbox:
    ds_osnap['sigma0'] = gsw.sigma0(SA=ds_osnap['sal'], CT=ds_osnap['temp'])

    # Define potential density bins:
    sigma0_bins = np.arange(20, 29, 0.01)

    # Compute Total OSNAP diapycnal overturning stream function:
    ds_osnap['moc_total'] = compute_section_moc_tracer(ds=ds_osnap,
                                                    tracer_name='sigma0',
                                                    tracer_bins=sigma0_bins,
                                                    dir='-1',
                                                    mask=None,
                                                    )

    # Determine station indexes for OSNAP East section:
    station_OWest_OEast = ds_osnap.station.where(ds_osnap.longitude <= -44).max()

    # OSNAP East diapycnal overturning stream function:
    mask_OEast = ds_osnap.station >= station_OWest_OEast
    ds_osnap['moc_east'] = compute_section_moc_tracer(ds=ds_osnap,
                                                        tracer_name='sigma0',
                                                        tracer_bins=sigma0_bins,
                                                        dir='-1',
                                                        mask=mask_OEast,
                                                        )

    # OSNAP West diapycnal overturning stream function:
    mask_OWest = ds_osnap.station < station_OWest_OEast
    ds_osnap['moc_west'] = compute_section_moc_tracer(ds=ds_osnap,
                                                        tracer_name='sigma0',
                                                        tracer_bins=sigma0_bins,
                                                        dir='-1',
                                                        mask=mask_OWest,
                                                        )
    print("Completed: Computed eORCA025 ERA-5 v1 OSNAP AMOC stream functions.")

    # === Save OSNAP Section to NetCDF === #
    ds_osnap.to_netcdf(out_file, unlimited_dims=['time'])
    print(f"Completed: Saved eORCA025 ERA-5 v1 OSNAP data to: {out_file}.")


if __name__ == "__main__":
    # Define paths to NEMO grid files:
    domain_path = "/dssgfs01/scratch/npd/simulations/Domains/eORCA025/domain_cfg.nc"
    t_paths= sorted(glob.glob("/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1/eORCA025_ERA5_1m_grid_T_197*.nc"))
    v_paths = sorted(glob.glob("/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1/eORCA025_ERA5_1m_grid_V_197*.nc"))
    u_paths = sorted(glob.glob("/dssgfs01/scratch/npd/simulations/eORCA025_ERA5_v1/eORCA025_ERA5_1m_grid_U_197*.nc"))

    T_paths = {'temp': t_paths, 'sal': t_paths}
    U_paths = {'uo': u_paths, 'uo_eiv': u_paths, 'e3u': u_paths}
    V_paths = {'vo': v_paths, 'vo_eiv': v_paths, 'e3v': v_paths}

    # Define output file path:
    out_file = "/dssgfs01/scratch/npd/diagnostics/eORCA025_ERA5_v1/osnap_ERA5_v1/eORCA025_ERA5_v1_osnap_1976-01_2024-12.nc"

    # Extract OSNAP section and compute MOC:
    extract_osnap_moc(domain_path=domain_path,
                      T_paths=T_paths,
                      U_paths=U_paths,
                      V_paths=V_paths,
                      out_file=out_file,
                      )
