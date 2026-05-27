"""
extend_ERA5_preprocessed_latest.py

Description: This script extends the valid_time dimension of the latest month of ERA5 preprocessed
data to ensure compatibility with fldread routine in NEMO SBC - which requires a full month of data.

Missing time steps are filled with NaN values.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
"""
# -- Import Dependencies -- #
import os
import logging
import xarray as xr
import numpy as np


# -- Define main function to extend ERA5 preprocessed data -- #
def main(filedir_latest: str = "/dssgfs01/scratch/npd/forcing/ERA5/preprocessed_latest",
         filedir_t2m_adj: str = "/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj",
         year: int = 2026,
         month: int = 5,
         variable_list: list = ["2m_temperature", "2m_dewpoint_temperature",
                                "mean_sea_level_pressure",
                                "mean_surface_downward_long_wave_radiation_flux",
                                "mean_surface_downward_short_wave_radiation_flux",
                                "mean_snowfall_rate", "mean_total_precipitation_rate",
                                "10m_u_component_of_wind", "10m_v_component_of_wind"]
         ) -> None:
    """
    Extend the valid_time dimension of the latest month of ERA5 preprocessed data.

    Parameters:
    -----------
    filedir_latest: str
        Directory of latest preprocessed ERA5 data to be extended.
    filedir_t2m_adj: str
        Directory of adjusted ERA5 2m temperature data to be extended.
    year: int
        Year of the latest month of ERA5 data to be extended.
    month: int
        Month of the latest month of ERA5 data to be extended.
    variable_list: list
        List of ERA5 variables to extend.
    """
    # --- Configure Logging --- #
    logging.basicConfig(
        filename="extend_ERA5_preprocessed_latest.log",
        encoding="utf-8",
        filemode="w",
        format="{asctime} - {levelname} - {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M",
        level=logging.INFO,
        )
    
    # -- Define Variable Mapping -- #
    var_map = {
        "2m_temperature": "t2m",
        "2m_dewpoint_temperature": "d2m",
        "mean_sea_level_pressure": "msl",
        "mean_surface_downward_long_wave_radiation_flux": "msdwlwrf",
        "mean_surface_downward_short_wave_radiation_flux": "msdwswrf",
        "mean_snowfall_rate": "msr",
        "mean_total_precipitation_rate": "mtpr",
        "10m_u_component_of_wind": "u10",
        "10m_v_component_of_wind": "v10"
    }

    logging.info(f"=== Extending ERA5 preprocessed data for {year}-{month:02d} ===")
    for var in variable_list:
        logging.info(f"-> In Progress: Extending ERA5 preprocessed variable: {var}")
        # -- Define Target Directory -- #
        if var == "2m_temperature":
            filedir = f"{filedir_t2m_adj}/{year}"
            filedir_ref = f"{filedir_t2m_adj}/{year - 1}"
            time_name = "time"
        else:
            filedir = f"{filedir_latest}/{year}/{var}"
            filedir_ref = f"{filedir_latest.replace('preprocessed_latest', 'preprocessed')}/{year - 1}/{var}"
            time_name = "valid_time"

        # -- Check if Output Filepath Exists -- #
        filepath_out = f"{filedir}/{var}_{year}-{month:02d}_extended.nc"
        if os.path.isfile(filepath_out):
            logging.info(f"Skipping file: {filepath_out} already exists.")
            continue

        # -- Open Reference ERA5 Dataset -- #
        filepath_ref = f"{filedir_ref}/{var}_{year - 1}-{month:02d}.nc"
        ds_ref = xr.open_dataset(filepath_ref)
        logging.info(f"-> Completed: Opened reference ERA5 dataset: {filepath_ref}")

        # -- Open Partial ERA5 Dataset -- #
        filepath_partial = f"{filedir}/{var}_{year}-{month:02d}.nc"
        ds_partial = xr.open_dataset(filepath_partial)
        logging.info(f"-> Completed: Opened latest partial ERA5 dataset: {filepath_partial}")

        # -- Define Empty ERA5 Dataset for Missing Time Steps -- #
        if len(ds_partial[time_name]) == len(ds_ref[time_name]):
            raise ValueError(f"ERA5 preprocessed variable {var} for {year}-{month:02d} already contains a full month of data.")
        nvt = len(ds_partial[time_name])
        # Select outstanding time steps from reference dataset and fill with NaN values:
        ds_ref = ds_ref.isel({time_name: slice(nvt, None)})
        ds_ref[var_map[var]].data = np.full_like(ds_ref[var_map[var]].data, fill_value=np.nan)

        # -- Concatenate Partial and Reference Datasets -- #
        ds_out = xr.concat([ds_partial, ds_ref], dim=time_name)

        # -- Write Extended Dataset to File -- #
        ds_out.to_netcdf(filepath_out)
        logging.info(f"-> Completed: Saved extended ERA5 preprocessed variable {var} to netCDF file: {filepath_out}")
        logging.info("  ======  ")

if __name__ == "__main__":
    # ========== Input Arguments ========== #
    # Directory of latest preprocessed ERA5 data to be extended:
    filedir_latest = "/dssgfs01/scratch/npd/forcing/ERA5/preprocessed_latest"

    # Directory of adjusted ERA5 2m temperature data to be extended:
    filedir_t2m_adj = "/dssgfs01/scratch/npd/forcing/ERA5_t2m_adj"

    # Define year and month of the latest month of ERA5 data to be extended:
    year = 2026
    month = 5

    # Define list of ERA5 variables to extend:
    variable_list = ["2m_temperature", "2m_dewpoint_temperature",
                    "mean_sea_level_pressure",
                    "mean_surface_downward_long_wave_radiation_flux",
                    "mean_surface_downward_short_wave_radiation_flux",
                    "mean_snowfall_rate", "mean_total_precipitation_rate",
                    "10m_u_component_of_wind", "10m_v_component_of_wind"
                    ]
    # ========== Input Arguments ========== #
    
    # -- Run Main Function -- #
    main(filedir_latest=filedir_latest,
         filedir_t2m_adj=filedir_t2m_adj,
         year=year,
         month=month,
         variable_list=variable_list
         )