#!/usr/bin/env bash
#SBATCH --job-name=OSNAP_eORCA025_ERA5v1
#SBATCH --partition=compute
#SBATCH --time=03:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=1

# ------------------------------------------------------------------------------
# run_extract_osnap_eORCA025.sh
# 
# Description: Script to calculate OSNAP AMOC diagnostics for eORCA025 Near-
# Present Day (NPD) simulation.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
#
# Created On: 2025-05-20
# ------------------------------------------------------------------------------

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_nemo_cookbook

# -- Run METRIC -- #
echo "In Progress: Calculating eORCA025 ERA-5 v1 OSNAP AMOC Diagnostics..."

# Call OSNAP extraction script:
python3 extract_osnap_eORCA025_ERA5v1.py

echo "Completed: Calculated eORCA025 ERA-5 v1 OSNAP AMOC Diagnostics."
