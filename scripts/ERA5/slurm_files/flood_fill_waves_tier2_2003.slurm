#!/bin/bash
#SBATCH --job-name=flooding
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-socket=8

module load netCDF/4.8.1-iimpi-2021b
module load NCO/5.0.3-foss-2021b

# Location of directories:
export cdodir=/dssgfs01/working/acc/TOOLS/bin
#export maskfile='/dssgfs01/scratch/npd/forcing/ERA5/original/land_sea_mask_binary_0.nc'
export indir='/dssgfs01/scratch/npd/forcing/ERA5/original/'
export outdir='/dssgfs01/scratch/npd/forcing/ERA5/preprocessed/'

# Input files:
year=2003
infiles=(${indir}${year}/{significant_height_of_total_swell,significant_height_of_wind_waves,mean_wave_period_based_on_first_moment_for_swell,mean_wave_period_based_on_first_moment_for_wind_waves,mean_direction_of_total_swell,mean_direction_of_wind_waves,mean_period_of_total_swell,mean_period_of_wind_waves,}/*.nc)
inf=${#infiles[*]}

# Find files not yet computed:
compfiles=()
for ff in $( seq 0 $(($inf-1)) )
do
  nfile=${infiles[$ff]##*/}
  outfile=${outdir}${nfile: -10: -6}/${nfile::-11}/${nfile}
  # If file doesn't exist add to list to be computed:
  if [[ -f "$outfile" ]]
  then
    echo "$outfile exists."
  else
    mkdir -p ${outdir}${nfile: -10: -6}/${nfile::-11}
    compfiles+=($nfile)
  fi
done
cnf=${#compfiles[*]}
echo Number of files to process: $cnf
echo Infiles: ${compfiles[*]}

export floodtemp=${outdir}waveFLOODTEMP${year}
export masktemp=${outdir}waveMASKTEMP${year}
mkdir -p $floodtemp
mkdir -p $masktemp

xargs -n 1 -P 16 <<< ${compfiles[*]} sh -c 'ncks -4 ${indir}${1: -10: -6}/${1::-11}/${1} ${masktemp}/${1}' bash
xargs -n 1 -P 16 <<< ${compfiles[*]} sh -c '${cdodir}/cdo -k auto -z zip_1 fillmiss3 ${masktemp}/${1} ${floodtemp}/${1}' bash
xargs -n 1 -P 12 <<< ${compfiles[*]} sh -c 'ncks -4 --mk_rec_dmn time --cnk_dmn longitude,24 --cnk_dmn latitude,24 --cnk_dmn time,1 --dfl_lvl 1 ${floodtemp}/${1} ${outdir}${1: -10: -6}/${1::-11}/${1}' bash

rm -rf $floodtemp
rm -rf $masktemp
