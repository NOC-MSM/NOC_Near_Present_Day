#!/bin/bash
#SBATCH --job-name=preprocess_ERA5_original
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-socket=8

# ------------------------------------------------------------
# run_preprocess_ERA5_forcing_original.slurm
# Description: 
# Perform flood filling, chunking and variable renaming for
# ERA-5 1-hourly atmospheric forcing fields & save to
# preprocessed directory.
# 
# Created By: Jenny Mecking (jmecki@noc.ac.uk)
# Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# ------------------------------------------------------------
# -- Load Dependencies -- #
module load netCDF/4.9.2-iimpi-2023a
module load NCO/5.1.9-foss-2023a

# -- Define Directories -- #
export cdodir=/dssgfs01/working/atb299/TOOLS/bin/
export maskfile='/dssgfs01/scratch/npd/forcing/ERA5/original/land_sea_mask_binary_0.nc'
export indir='/dssgfs01/scratch/npd/forcing/ERA5/original/'
export outdir='/dssgfs01/scratch/npd/forcing/ERA5/preprocessed/'

# -- Define Input Files -- #
# TO DO: Update from start_year below to present.
year=2025
infiles=(${indir}${year}/{10m_u_component_of_wind,10m_v_component_of_wind,2m_dewpoint_temperature,2m_temperature,mean_sea_level_pressure,mean_snowfall_rate,mean_surface_downward_short_wave_radiation_flux,mean_surface_downward_long_wave_radiation_flux,mean_total_precipitation_rate}/*.nc)
inf=${#infiles[*]}

# -- Find Unprocessed Files -- #
compfiles=()
for ff in $( seq 0 $(($inf-1)) )
do
  nfile=${infiles[$ff]##*/}
  outfile=${outdir}${nfile: -10: -6}/${nfile::-11}/${nfile}
  # If file doesn't exist add to list to be computed:
  if [[ -f "$outfile" ]]
  then
    echo "Skipping file: $outfile already exists."
  else
    mkdir -p ${outdir}${nfile: -10: -6}/${nfile::-11}
    compfiles+=($nfile)
  fi
done
cnf=${#compfiles[*]}
echo Number of Files to Flood Fill: $cnf
echo File List: ${compfiles[*]}

# -- Create Temporary Directories -- #
export masktemp=${outdir}MASKTEMP${year}
export floodtemp=${outdir}FLOODTEMP${year}
mkdir -p $masktemp
mkdir -p $floodtemp

# -- Flood Fill ERA-5 Forcing Fields -- #
echo "================================================="
echo "In Progress: Flood filling ERA-5 forcing files..."
xargs -n 1 -P 16 <<< ${compfiles[*]} sh -c '${cdodir}/cdo -k auto -z zip_1 ifthen $maskfile ${indir}${1: -10: -6}/${1::-11}/${1} ${masktemp}/${1}' bash
xargs -n 1 -P 16 <<< ${compfiles[*]} sh -c '${cdodir}/cdo -k auto -z zip_1 fillmiss3 ${masktemp}/${1} ${floodtemp}/${1}' bash
echo "Completed: Flood filled ERA-5 forcing files."

echo "================================================="
echo "In Progress: Rechunking ERA-5 forcing files..."
xargs -n 1 -P 6 <<< ${compfiles[*]} sh -c 'ncks --mk_rec_dmn time --cnk_dmn longitude,24 --cnk_dmn latitude,24 --cnk_dmn time,1 --dfl_lvl 1 ${floodtemp}/${1} ${outdir}${1: -10: -6}/${1::-11}/${1}' bash
echo "Completed: Rechunked ERA-5 forcing files."

# -- Renaming Variables ERA-5 Forcing Fields -- #
if [ $year -ge 2025 ]; then
  echo "================================================="
  echo "In Progress: Renaming surface flux variables in $year ERA-5 forcing files..."
    for ff in $( seq 0 $(($inf-1)) )
    do
      nfile=${compfiles[$ff]}
      npath=${outdir}${nfile: -10: -6}/${nfile::-11}/${nfile}
      case $nfile in
        mean_snowfall_rate*)
          ncrename -v avg_tsrwe,msr ${npath}
          ;;
        mean_surface_downward_short_wave_radiation_flux*)
          ncrename -v avg_sdswrf,msdwswrf ${npath}
          ;;
        mean_surface_downward_long_wave_radiation_flux*)
          ncrename -v avg_sdlwrf,msdwlwrf ${npath}
          ;;
        mean_total_precipitation_rate*)
          ncrename -v avg_tprate,mtpr ${npath}
          ;;
        *)
          continue
          ;;
      esac
    done
  echo "Completed: Renamed Variables in ERA-5 forcing files."
fi

# -- Remove Temporary Directories -- #
rm -rf $masktemp
rm -rf $floodtemp
