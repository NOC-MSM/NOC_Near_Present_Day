#!/bin/bash
#SBATCH --job-name=NPD_12
#SBATCH --account=n01-CLASS
#SBATCH --time=17:10:00
#SBATCH --nodes=134
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:cyclic

set -u
set -e
set -o pipefail


module load cray-mpich/8.1.23 
#module swap craype-network-ofi craype-network-ucx
#module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1 
module load cray-netcdf-hdf5parallel/4.9.0.1

export NEMO_EXP_DIR="/work/n01/n01/acc/NEMO/NPD/NOC_Near_Present_Day/nemo/cfgs/GLOBAL_QCO/eORCA12"
export XIO_HOME="/work/n01/shared/nemo/XIOS3_trunk_Cray"
export OMP_NUM_THREADS="1"
export I_MPI_SHM="off"
export XIOS_PROCESSES="192"
export NEMO_PROCESSES="5504"
#export UCX_UD_MLX5_TIMEOUT="2m"

cd "$NEMO_EXP_DIR"

if [ -f namelist_cfg ]
then
 nn_itend=`grep nn_itend namelist_cfg | grep -v '!' | sed -e 's/.*=\s*//' -e 's/,//'`
else
 nn_itend=1
fi

# Link nemo
sbcast --compress=none "$NEMO_EXP_DIR"/../BLD/bin/nemo.exe /tmp/nemo
sbcast --compress=none "$XIO_HOME"/bin/xios_server.exe /tmp/xios_server.exe
# Run NEMO

echo "Commence srun.."

export SBATCH_NO_KILL=1
srun --mem-bind=local \
--het-group=0 --nodes=48 --ntasks=192 --ntasks-per-core=1 --ntasks-per-node=4 --cpus-per-task=16 /tmp/xios_server.exe \
: --het-group=1 --nodes=86 --ntasks=5504 --ntasks-per-core=1 --ntasks-per-node=64 --cpus-per-task=2  /tmp/nemo &
srun_pid=$!
./time_step_chker $srun_pid $nn_itend testrun.slurm &
   jobs
   wait %1

      echo -n "Rerun required? "
      if [ -f failed_to_start ] ; then  echo yes ; else echo no ; fi
      if [ -f failed_to_start ]
      then
       echo "Go again"
       kill %2
       rm failed_to_start
       touch time.step
       srun --mem-bind=local \
       --het-group=0 --nodes=48 --ntasks=192 --ntasks-per-core=1 --ntasks-per-node=4 --cpus-per-task=16 /tmp/xios_server.exe \
       : --het-group=1 --nodes=86 --ntasks=5504 --ntasks-per-core=1 --ntasks-per-node=64 --cpus-per-task=2  /tmp/nemo &
       srun_pid=$!
       ./time_step_chker $srun_pid $nn_itend testrun.slurm &
       jobs
       wait %3
      fi
