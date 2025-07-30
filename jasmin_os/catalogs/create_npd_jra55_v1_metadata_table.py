"""
create_npd_jra55v1_metadata_table.py

Description: This script creates a metadata table for the National
Oceanography Centre's Near-Present-Day (NPD) JRA55-do version 1
outputs. The metadata table is saved in a single CSV file.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
Creation Date: 2025-03-03
"""
# -- Import required packages -- #
import glob
import pandas as pd
import xarray as xr

# -- Create catalog DataFrame -- #
# Define empty dataframe:
df_inv = pd.DataFrame(columns=['variable', 'standard_name', 'long_name', 'units', 'dims', 'model', 'grid', 'freq', 'url'])

# Iterate over models, grids, and output frequencies:
for model in ['eORCA1', 'eORCA025']:
    for grid in ['T', 'U', 'V', 'W']:
        for freq in ['5d', '1m', '1y']:
            # Ignore eORCA1 5d frequency since this is not available:
            if (model == 'eORCA1') & (freq == '5d'):
                continue
            else:
                # Open example dataset to get variable metadata:
                filepaths = glob.glob(f'/dssgfs01/scratch/otooth/npd_data/simulations/{model}_JRA55/exp_npd_v1/{model}_{freq}_grid_{grid}*.nc')
                filepaths.sort()
                ds = xr.open_dataset(filepaths[0])
                # Iterate over available variables:
                for var in ds.data_vars:
                    attrs = ds[var].attrs
                    # Ignore variables without standard_name attribute (e.g. coordinate dimensions):
                    if 'standard_name' not in attrs:
                        continue
                    else:
                        # Create JASMIN OS read-only URL:
                        url = f'https://noc-msm-o.s3-ext.jc.rl.ac.uk/npd-{model.lower()}-jra55v1/{grid}{freq}/{var}'
                        # Append variable metadata to catalog dataframe:
                        df_inv.loc[len(df_inv)] = [var, attrs['standard_name'].lower(), attrs['long_name'].lower(), attrs['units'], ds[var].dims, model, grid, freq, url]

# Save catalog DataFrame to CSV file:
df_inv.to_csv('.../NOC_Near_Present_Day/jasmin_os/intake/npd_jra55_v1_catalog.csv', index=False)