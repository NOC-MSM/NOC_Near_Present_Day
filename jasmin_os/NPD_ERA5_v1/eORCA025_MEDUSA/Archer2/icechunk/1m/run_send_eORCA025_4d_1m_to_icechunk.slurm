#!/bin/bash
#SBATCH --job-name=xfer_4d_1m_eORCA025
#SBATCH --time=12:00:00
#SBATCH --ntasks=22
#SBATCH --mem=65GB

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# ==============================================================
# run_send_eORCA025_4d_1m_to_icechunk.slurm
#
# Description: SLURM script to send the eORCA025-MEDUSA ERA5v1
# monthly mean 4-D variables to icechunk repositories.
#
# Created By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# Created On: 2026-03-23
#
# ==============================================================
set -euo pipefail

# -- Input arguments to OceanDataStore -- #
# Filepath to eORCA025 ancillary file:
filepath_grid=/work/n01/n01/julpal/NEMO/NPD_cray_newer/nemo/cfgs/GLOBAL_QCO_MEDUSA/eORCA025/INPUT/domcfg.nc

# Filepath to eORCA025 monthly mean output files:
filedir=/work/n01/n01/julpal/NEMO/NPD_cray_newer/nemo/cfgs/GLOBAL_QCO_MEDUSA/eORCA025/OUTPUT
# NEMO + SI3 Output:
filepath_gridT=$filedir/eORCA025_1m_grid_T_*.nc
filepath_gridU=$filedir/eORCA025_1m_grid_U_*.nc
filepath_gridV=$filedir/eORCA025_1m_grid_V_*.nc
filepath_gridW=$filedir/eORCA025_1m_grid_W_*.nc
# MEDUSA Biogeochemical Output:
filepath_gridP=$filedir/eORCA025_1m_ptrc_T_*.nc
filepath_gridD=$filedir/eORCA025_1m_diad_T_*.nc

# Filepath to JASMIN OS credentials:
store_credentials_json=./credentials/jasmin_os_credentials.json

# Filepath to Dask configuration:
dask_config_json=./dask_config.json

# Define bucket:
bucket=npd-eorca025-medusa-era5v1

# -- Python Environment -- #
# Activate miniconda environment:
source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate env_jasmin_os

# === Send eORCA025-MEDUSA ERA5v1 4-D T1m variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA ERA-5 v1 4-D T1m variables to Icechunk..."

# Define 4-dimensional T1m variables:
variables="time_centered_bounds time_counter_bounds e3t thetao_con so_abs rsdo rsdoabsorb"
# variables="time_centered_bounds time_counter_bounds e3t thetao_con so_abs rsdo vohflisf vohfcisf vowflisf ocontempadvect ocontemppmdiff ocontempdiff ttrd_evd_li rsdoabsorb ttrd_bbl_li ttrd_atf_li ocontemptend osaltadvect osaltpmdiff osaltdiff strd_evd_li osalttend strd_bbl_li strd_atf_li"

ods send_to_icechunk -f $filepath_gridT -c $store_credentials_json -b $bucket -p T1m_4d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":1440, "y":1206}' -dc $dask_config_json -v $variables \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 T1m (1976-01 - 2025-07)" -vc || { echo "Error: ods send_to_icechunk failed."; exit 1; }

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 4-D T1m monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 4-D U1m variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA ERA-5 v1 4-D U1m variables to Icechunk..."

# Define 4-dimensional U1m variables:
variables="time_centered_bounds time_counter_bounds e3u uo uo_eiv u2o umo"

ods send_to_icechunk -f $filepath_gridU -c $store_credentials_json -b $bucket -p U1m_4d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamu", "nav_lat":"gphiu"}' \
                     -cs '{"x":1440, "y":1206}' -dc $dask_config_json -v $variables \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 U1m (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 4-D U1m monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 4-D V1m variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA ERA-5 v1 4-D V1m variables to Icechunk..."

# Define 4-dimensional variables:
variables="time_centered_bounds time_counter_bounds e3v vo vo_eiv v2o vmo"

ods send_to_icechunk -f $filepath_gridV -c $store_credentials_json -b $bucket -p V1m_4d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamv", "nav_lat":"gphiv"}' \
                     -cs '{"x":1440, "y":1206}' -dc $dask_config_json -v $variables \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 V1m (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 4-D V1m monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 4-D W1m variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA ERA-5 v1 4-D W1m variables to Icechunk..."

# Define 4-dimensional variables:
variables="time_centered_bounds time_counter_bounds e3w wmo wo w2o obvfsq avt_evd difvho difvso difvmo diftrto"

ods send_to_icechunk -f $filepath_gridW -c $store_credentials_json -b $bucket -p W1m_4d \
                     -gf $filepath_domain -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":720, "y":603, "depthw":25}' -dc $dask_config_json -v $variables \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 W1m (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 4-D W1m monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 4-D P1m variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA 4-D P1m (i.e., ptrc_T) variables to Icechunk..."

# Define 4-dimensional P1m variables:
variables_gridP="time_centered_bounds time_counter_bounds CHN CHD PHN PHD ZMI ZME DIN SIL FER DET PDS DTC DIC ALK OXY"

ods send_to_icechunk -f $filepath_gridP -c $store_credentials_json -b $bucket -p P1m_4d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":720, "y":603, "deptht":25}' -v $variables_gridP  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 P1m (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 4-D P1m monthly mean fields to Icechunk."

# === Send eORCA025-MEDUSA ERA5v1 4-D D1m variables to Icechunk === #
echo "In Progress: Sending eORCA025-MEDUSA 4-D D1m (i.e., diad) variables to Icechunk..."

# Define 4-dimensional D1m variables:
# Subset of diagnostic Variables:
variables_gridD="time_centered_bounds time_counter_bounds INVTN INVTSI INVTFE PRN MPN PRD MPD OPAL OPALDISS GMIPN GMID MZMI GMEPN GMEPD GMEZMI GMED MZME DETN MDET AEOLIAN BENTHIC SCAVENGE PN_JLIM PN_NLIM PN_FELIM PD_JLIM PD_NLIM PD_FELIM PD_SILIM PDSILIM2 SDT__100 SDT__200 SDT__500 SDT_1000 TOTREG_N TOTRG_SI REG__100 REG__200 REG__500 REG_1000 FASTN FASTSI FASTFE FASTC FASTCA FDT__100 FDT__200 FDT__500 FDT_1000 RG__100F RG__200F RG__500F RG_1000F FDS__100 FDS__200 FDS__500 FDS_1000 RGS_100F RGS_200F RGS_500F RGS1000F REMINN REMINSI REMINFE REMINC REMINCA SEAFLRN SEAFLRSI SEAFLRFE SEAFLRC SEAFLRCA MED_QSR MED_XPAR INTFLX_N INTFLX_SI INTFLX_FE TOT_ISEDN TOT_ISEDSI TOT_ISEDFE TOT_ISEDC TOT_ISEDCA TOT_OSEDN TOT_OSEDSI TOT_OSEDFE TOT_OSEDC TOT_OSEDCA INT_PN INT_PD ML_PRN ML_PRD OCAL_CCD OCAL_LVL FE_0000 FE_0100 FE_0200 FE_0500 FE_1000 MED_XZE WIND ATM_PCO2 OCN_PH OCN_PCO2 OCNH2CO3 OCN_HCO3 OCN_CO3 CO2FLUX FGCO2 OM_CAL OM_ARG TCO2 TALK KW660 ATM_PP0 O2FLUX O2SAT CAL_CCD ARG_CCD SFR_OCAL SFR_OARG N_PROD N_CONS C_PROD C_CONS O2_PROD O2_CONS O2_ANOX RR_0100 RR_0500 RR_1000 IBEN_N IBEN_FE IBEN_C IBEN_SI IBEN_CA OBEN_N OBEN_FE OBEN_C OBEN_SI OBEN_CA BEN_N BEN_FE BEN_C BEN_SI BEN_CA DETC SDC__100 SDC__200 SDC__500 SDC_1000 INVTC INVTALK INVTO2 LYSO_CA COM_RESP PN_LLOSS PD_LLOSS ZI_LLOSS ZE_LLOSS ZI_MES_N ZI_MES_D ZI_MES_C ZI_MESDC ZI_EXCR ZI_RESP ZI_GROW ZE_MES_N ZE_MES_D ZE_MES_C ZE_MESDC ZE_EXCR ZE_RESP ZE_GROW MDETC GMIDC GMEDC INT_ZMI INT_ZME INT_DET INT_DTC ATM_XCO2 OCN_FCO2 ATM_FCO2 OCN_RHOSW OCN_SCHCO2 OCN_KWCO2 OCN_K0 CO2STARAIR OCN_DPCO2 CHL_MLD RUNOFF RIV_N RIV_SI RIV_C RIV_ALK TPP3 DETFLUX3 REMIN3N PH3 OM_CAL3 MED_XPAR3D O2_PROD3 O2_CONS3 BDDTO2 epC100 epCALC100 epN100 epSI100 INTDISSIC INTDISSIN INTDISSISI INTTALK O2min ZO2min"

ods send_to_icechunk -f $filepath_gridD -c $store_credentials_json -b $bucket -p D1m_4d \
                     -gf $filepath_grid -uc '{"nav_lon":"glamt", "nav_lat":"gphit"}' \
                     -cs '{"x":720, "y":603, "deptht":25}' -v $variables_gridD  -dc $dask_config_json \
                     -br "main" -cm "Added eORCA025-MEDUSA ERA5v1 D1m (1976-01 - 2025-07)" -vc

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 4-D D1m monthly mean fields to Icechunk."

echo "Completed: Sent eORCA025-MEDUSA ERA-5 v1 4-D monthly mean fields to Icechunk."