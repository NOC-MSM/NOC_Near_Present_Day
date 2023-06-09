#!/bin/bash
#SBATCH --job-name=rebuild
#SBATCH --time=02:20:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks-per-core=1
##SBATCH --ntasks-per-node=16
##SBATCH --ntasks-per-socket=8
#SBATCH --ntasks-per-node=32
#SBATCH --ntasks-per-socket=16

module load NEMO/prg-env

#############################################
# Run rebuild_nemo on all sets of output files in the OUTPUT directory
# - Uses xargs to run N instances of rebuild_nemo in parallel "-P 16" 
# - Note rebuild_nemo is also parallel "-p 4"
#
# A. Blaker 12/05/2021
#############################################

export I_MPI_SHM=icx
echo "I_MPI_SHM = " $I_MPI_SHM
export OMP_NUM_THREADS=4

#============================
TIME1=`date +%s`

# This is dangerous - do not run whilst simulation is in progress
ff=`ls ../OUTPUT/*scalar*000.nc`
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 73 <<<$ff sh -c 'mv ${1} ${1%_*}.nc; rm ${1%_*}_????.nc' bash
fi

ff=`ls ../OUTPUT/*_1y*000.nc`
##     ^^^^^^^^^^^^^^^^^^^^^
##     Change this to point to your output folder/location
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 8 <<<$ff sh -c '/dssgfs01/working/atb299/REBUILD_NEMO/rebuild_nemo -p 4 -d 1 ${1%_*} 24 &> `basename ${1%_*}`.rpt' bash
#                                                                                                     ^^
#                                                                                Set this to match the # XIOS servers you are using
fi

ff=`ls ../OUTPUT/*_1m*000.nc`
##     ^^^^^^^^^^^^^^^^^^^^^
##     Change this to point to your output folder/location
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 8 <<<$ff sh -c '/dssgfs01/working/atb299/REBUILD_NEMO/rebuild_nemo -p 4 -d 1 ${1%_*} 24 &> `basename ${1%_*}`.rpt' bash
#                                                                                                     ^^
#                                                                                Set this to match the # XIOS servers you are using
fi

ff=`ls ../OUTPUT/*_1d*000.nc`
##     ^^^^^^^^^^^^^^^^^^^^^
##     Change this to point to your output folder/location
if [ ! -z "$ff" ]; then
  xargs -n 1 -P 8 <<<$ff sh -c '/dssgfs01/working/atb299/REBUILD_NEMO/rebuild_nemo -p 4 -d 1 ${1%_*} 24 &> `basename ${1%_*}`.rpt' bash
#                                                                                                    ^^
#                                                                                Set this to match the # XIOS servers you are using
fi

#============================
# Job timing

TIME2=`date +%s`
DIFFSEC=`expr ${TIME2} - ${TIME1}`
echo Took ${DIFFSEC} seconds.
echo Took `date +%H:%M:%S -ud @${DIFFSEC}`

