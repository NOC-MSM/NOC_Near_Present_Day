#!/bin/bash
#SBATCH --job-name=rebuild
#SBATCH --time=02:20:00
#SBATCH --account=n01-ORCHESTRA
#SBATCH --partition=standard
#SBATCH --qos=standard
##SBATCH --qos=short

#SBATCH --nodes=1

module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7

#############################################
# Run rebuild_nemo on all sets of output files in the OUTPUTS directory
# - Uses xargs to run N instances of rebuild_nemo in parallel "-P 16" 
# - Note rebuild_nemo is also parallel "-p 4"
#
# A. Blaker 12/05/2021
#############################################

rebuild_nemo='/work/n01/n01/atb299/TOOLS/tools/REBUILD_NEMO/BLD/bin/rebuild_nemo'

export OMP_NUM_THREADS=1
#============================
TIME1=`date +%s`

ff=`ls ../OUTPUT/*scalar*000.nc`
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 73 <<<$ff sh -c 'mv ${1} ${1%_*}.nc; rm ${1%_*}_????.nc' bash
fi

ff=`ls ../OUTPUT/*_1d*000.nc`
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 32 <<<$ff sh -c '/work/n01/n01/atb299/TOOLS/tools/REBUILD_NEMO/rebuild_nemo -p 4 -d 1 ${1%_*} 48 &> `basename ${1%_*}`.rpt' bash
fi

ff=`ls ../OUTPUT/*_1m*000.nc`
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 16 <<<$ff sh -c '/work/n01/n01/atb299/TOOLS/tools/REBUILD_NEMO/rebuild_nemo -p 4 -d 1 ${1%_*} 48 &> `basename ${1%_*}`.rpt' bash
fi

ff=`ls ../OUTPUT/*_5d*000.nc`
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 8 <<<$ff sh -c '/work/n01/n01/atb299/TOOLS/tools/REBUILD_NEMO/rebuild_nemo -p 4 -d 1 ${1%_*} 48 &> `basename ${1%_*}`.rpt' bash
fi

#============================
# Job timing

TIME2=`date +%s`
DIFFSEC=`expr ${TIME2} - ${TIME1}`
echo Took ${DIFFSEC} seconds.
echo Took `date +%H:%M:%S -ud @${DIFFSEC}`

