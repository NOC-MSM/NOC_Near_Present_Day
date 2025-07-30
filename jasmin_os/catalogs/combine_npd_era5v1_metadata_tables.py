"""
combine_npd_era5v1_metadata_tables.py

Description: This script combines the metadata tables for the National
Oceanography Centre's Near-Present-Day (NPD) ERA-5 version 1
outputs. The final table is saved as a single CSV file.

Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
Creation Date: 2025-07-30
Last Modified: 2025-07-30
"""
# -- Import required packages -- #
import pandas as pd

# -- Define the path to the metadata files -- #
eorca1_path = "./npd_eorca1_era5v1_catalog.csv"
eorca025_path = "./npd_eorca025_era5v1_catalog.csv"
eorca12_path = "./npd_eorca12_era5v1_catalog.csv"

# -- Load the metadata tables -- #
df_eorca1 = pd.read_csv(eorca1_path)
df_eorca025 = pd.read_csv(eorca025_path)
df_eorca12 = pd.read_csv(eorca12_path)

# -- Combine the metadata tables -- #
df_combined = pd.concat([df_eorca1, df_eorca025, df_eorca12], ignore_index=True)

# -- Store the combined metadata table as a CSV file -- #
output_path = "./npd_era5_v1_catalog.csv"
df_combined.to_csv(output_path, index=False)
print(f"Combined metadata table saved to: {output_path}")