"""
create_npd_era5v1_metadata_tables.py

Description: This script creates metadata tables for the National
Oceanography Centre's Near-Present-Day (NPD) ERA-5 version 1
outputs. The tables are saved as separate CSV files.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
Creation Date: 2025-03-03
Last Modified: 2025-07-30
"""
# -- Import required packages -- #
import glob
import pandas as pd
import xarray as xr

arch = 'Archer2'  # Set architecture type: 'Anemone' or 'Archer2'

# -- Anemone Metadata DataFrames -- #
if arch == 'Anemone':
    # Iterate over models, grids, and output frequencies:
    for model in ['eORCA1', 'eORCA025']:
        # Define empty dataframe:
        df_inv = pd.DataFrame(columns=['variable', 'standard_name', 'long_name', 'units', 'dims', 'model', 'grid', 'freq', 'url'])
        # Iterate over grid types and output frequencies:
        for grid in ['T', 'U', 'V', 'W']:
            for freq in ['5d', '1m', '1y']:
                # Ignore eORCA1 5d frequency since this is not available:
                if (model == 'eORCA1') & (freq == '5d'):
                    continue
                else:
                    # Open example dataset to get variable metadata:
                    filepaths = glob.glob(f'/dssgfs01/scratch/npd/simulations/{model}_ERA5_v1/{model}_ERA5_{freq}_grid_{grid}*.nc')
                    filepaths.sort()
                    ds = xr.open_dataset(filepaths[0])
                    # Iterate over available variables:
                    for var in ds.data_vars:
                        attrs = ds[var].attrs
                        # Ignore variables without standard_name attribute (e.g. coordinate dimensions):
                        if 'standard_name' not in attrs:
                            continue
                        else:
                            if model == 'eORCA1':
                                # Create JASMIN OS path for eORCA1 variables:
                                url = f'https://noc-msm-o.s3-ext.jc.rl.ac.uk/npd-{model.lower()}-era5v1/{grid}{freq}'
                            elif ds[var].ndim == 3:
                                # Create JASMIN OS path for 3D variables:
                                url = f'https://noc-msm-o.s3-ext.jc.rl.ac.uk/npd-{model.lower()}-era5v1/{grid}{freq}_3d'
                            elif ds[var].ndim == 4:
                                # Create JASMIN OS path for 4D variables:
                                url = f'https://noc-msm-o.s3-ext.jc.rl.ac.uk/npd-{model.lower()}-era5v1/{grid}{freq}_4d'

                            # Append variable metadata to catalog dataframe:
                            df_inv.loc[len(df_inv)] = [var, attrs['standard_name'].lower(), attrs['long_name'].lower(), attrs['units'], ds[var].dims, model, grid, freq, url]

        # Save model catalog DataFrame to CSV file:
        df_inv.to_csv(f'/dssgfs01/scratch/otooth/NOC_NPD/20250512_eORCA1_ERA5v1/NOC_Near_Present_Day/jasmin_os/catalogs/npd_{model.lower()}_era5v1_catalog.csv', index=False)

# -- Archer2 Metadata DataFrame -- #
elif arch == 'Archer2':
    # eORCA12 NPD model:
    model = 'eORCA12'
    # Define empty dataframe:
    df_inv = pd.DataFrame(columns=['variable', 'standard_name', 'long_name', 'units', 'dims', 'model', 'grid', 'freq', 'url'])
    # Iterate over grid types and output frequencies:
    for grid in ['T', 'U', 'V', 'W']:
        for freq in ['5d', '1m', '1y']:
            # Open example dataset to get variable metadata:
            filepaths = glob.glob(f'/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/{model}/TIDY/????/{model}_{freq}_grid_{grid}*.nc')
            filepaths.sort()
            ds = xr.open_dataset(filepaths[0])
            # Iterate over available variables:
            for var in ds.data_vars:
                attrs = ds[var].attrs
                # Ignore variables without standard_name attribute (e.g. coordinate dimensions):
                if 'standard_name' not in attrs:
                    continue
                else:
                    if ds[var].ndim == 3:
                        # Create JASMIN OS path for 3D variables:
                        url = f'https://noc-msm-o.s3-ext.jc.rl.ac.uk/npd-{model.lower()}-era5v1/{grid}{freq}_3d'
                    elif ds[var].ndim == 4:
                        # Create JASMIN OS path for 4D variables:
                        url = f'https://noc-msm-o.s3-ext.jc.rl.ac.uk/npd-{model.lower()}-era5v1/{grid}{freq}_4d'

                    # Append variable metadata to catalog dataframe:
                    df_inv.loc[len(df_inv)] = [var, attrs['standard_name'].lower(), attrs['long_name'].lower(), attrs['units'], ds[var].dims, model, grid, freq, url]

    # Save model catalog DataFrame to CSV file:
    df_inv.to_csv(f'/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/jasmin_os/catalogs/npd_{model.lower()}_era5v1_catalog.csv', index=False)
