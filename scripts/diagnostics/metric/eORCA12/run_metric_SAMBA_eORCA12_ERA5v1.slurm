#!/bin/bash
#SBATCH --job-name=metric_SAMBA_M1m_eORCA12
#SBATCH --time=02:00:00
#SBATCH --ntasks=26
#SBATCH --mem=52GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

GLOBIGNORE="*"
# ------------------------------------------------------------------------------
# run_metric_SAMBA_eORCA12.sh
# 
# Description: Script to run METRIC SAMBA 34.5S AMOC diagnostic for eORCA12
# Near-Present Day (NPD) simulation.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
#
# Created On: 2025-08-15
# ------------------------------------------------------------------------------

# -- Input arguments to ./run_metric -- #
# Filepaths to config file:
fpath_config=/work/n01/n01/otooth/NOC-NPD/NEMO_cfgs/20250312_eORCA12/NOC_Near_Present_Day/scripts/diagnostics/metric/eORCA12/config_SAMBA_eORCA12.ini

# Filepaths to eORCA12 monthly mean output files:
fdir="/work/n01/n01/atb299/NOC_NPD/20250228/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12/TIDY"

# -- Python Environment -- #
# Activate miniconda environment:
source /home2/home/n01/n01/otooth/miniforge3/bin/activate
conda activate env_npd_diags

# -- Run METRIC -- #
for yr in 198 199 200 202
do
    echo "In Progress: Calculating ${yr}0s SAMBA 34.5S AMOC Diagnostics..."
    fpath_T="${fdir}/${yr}*/eORCA12_1m_grid_T_*.nc"
    fpath_V="${fdir}/${yr}*/eORCA12_1m_grid_V_*.nc"
    fpath_U="${fdir}/${yr}*/eORCA12_1m_grid_U_*.nc"

    # Run METRIC .py script:
    ./run_metric.py -c $fpath_config -t $fpath_T -s $fpath_T -v $fpath_V -ssh $fpath_T -taux $fpath_U

    echo "Completed: Calculated ${yr}0s SAMBA 34.5S AMOC Diagnostics."

done
