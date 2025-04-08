#!/bin/bash
#SBATCH --job-name=NPD_025
#SBATCH --account=n01-CLASS
#SBATCH --time=00:20:00
#SBATCH --nodes=47
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --exclusive
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:cyclic


# ========================================================
# PARAMETERS TO SET
# time units used here for restart frequency and simulaion length
  TIME_UNITS=0 # 0=years ; 1=days ; 2=hours
# Restart/resubmission frequency (in TIME_UNITS)
  FREQRST=1
# job-step initial time step (0: infer from time.step)
# IT000 != 0 -> auto-resubmission is switched OFF
  IT000=1
#
# Simulation original starting time step (unchanged for LENGTHxTIME_UNITS)
  ITBEGIN=1
# Simulation length (in TIME_UNITS)
  LENGTH=1
# Name of this script (to resubmit)
  SCRIPTNAME=run_het_nemo2998_96x_47n_a2.slurm
#
# Use of nn_stocklist to trigger end-of-year restarts regardless of leap- or non-leap-
# years, means job-steps have to be limited to no more than 10 years. Force this if
# yearly time units are used:
#
if [ "$TIME_UNITS" -eq 0 ] ; then
    if [ "$FREQRST" -gt 10 ] ; then FREQRST=10 ; fi
fi
# ========================================================
## modules load
module load cray-mpich/8.1.23
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1
module load cray-netcdf-hdf5parallel/4.9.0.1
module load libfabric

export NEMO_EXP_DIR=$PWD
echo $NEMO_EXP_DIR
export OMP_NUM_THREADS="1"
export I_MPI_SHM="off"
unset I_MPI_PMI_LIBRARY
export FI_CXI_OPTIMIZED_MRS=false
# UD transport layer is "unreliable", so disable it by setting...
export UCX_TLS=rc,dc,self,sm

cd "$NEMO_EXP_DIR"

# Link nemo and xios -> Check consistency with arch file
ln -s /work/n01/shared/nemo/XIOS3_r2700/bin/xios_server.exe .
ln -s ../BLD/bin/nemo.exe .

# Broadcast executables
sbcast --compress=none "$NEMO_EXP_DIR"/nemo.exe /tmp/nemo
sbcast --compress=none "$NEMO_EXP_DIR"/xios_server.exe /tmp/xios_server.exe

# ========================================================
# Define useful functions
get_var () {
    VARNAME=$1
    NAMNAME=$2
    VAR=$(grep "^\s*${VARNAME}\s*=.*" "$NAMNAME" | sed "s/[', ]//g")
    VAR=${VAR%%!*}
    VAR=${VAR#*=}
    echo "$VAR"
}

edit_var () {
    VARNAME=$1
    VAR=$2
    NAMNAMES=${*:3}
    CMD=$(echo sed -i "\"s|^\s*${VARNAME}\s*=.*|\t${VARNAME} = ${VAR},|g\"" "$NAMNAMES")
    eval "$CMD"
}

get_timestep () {
    TIMENAME=$1
    TIMESTEP=$(grep -Eo '[0-9\.]+' "$TIMENAME")
    echo "$TIMESTEP"
}
date_from_restart() {
    # Retrieves the date variable from a restart file and also returns the appropriate year
    # for the itend calculation. If the restart is from a date in December, then the next year
    # is returned
    if [ -f $1 ] ; then
     DATEN=$(ncdump -v ndastp $1 | grep "ndastp =" | awk '{print $3}' )
     YRST=${DATEN:0:4}
     MTHN=${DATEN:4:2}
     DAY=${DATEN:6:2}
     if [ "$MTHN" == "12" ] ; then
      YNXT=$(( $YRST + 1 ))
     else
      YNXT=$YRST
      # YRST unchanged but restarting before the last month risks itend calculation not providing an end-of-year timestep
     fi
     echo $YRST $MTHN $DAY $YNXT
    else
     echo "1800 01 01 -1"
    fi
}
greg_length() {
    # Returns length of given year in days (Gregorian calendar)
    YR=$1
    if [ $(( $YR % 400 )) -eq 0 ] ; then
      echo 366
    elif [ $(( $YR % 100 )) -eq 0 ] ; then
      echo 365
    elif [ $(( $YR % 4 )) -eq 0 ] ; then
      echo 366
    else
     echo 365
    fi
}


# Get it000
if [ "$IT000" -eq 0 ]; then
    RESUBMIT=true
    if test -f "time.step"; then
        IT000=$(get_timestep time.step)
    fi
    IT000=$(( IT000 + 1 ))
else
    RESUBMIT=false
fi

# Get dt
DT=$(get_var rn_Dt namelist_cfg)
DT=${DT%.*}

# Edit namelists
if [ "$IT000" -gt 1 ]; then
    # Use restarts
    edit_var ln_rstart ".true." ./namelist_cfg
    edit_var ln_tsd_init ".false." ./namelist_cfg

    # Set start
    TSTEP=$((IT000 - 1))
    IT000=$((TSTEP + 1))
    edit_var "nn_it000" "$IT000" ./namelist_cfg

    # Set restart prefixes
    EXP=$(get_var cn_exp ./namelist_cfg)
    RST_DIR=$(get_var cn_ocerst_indir ./namelist_cfg)
    OCERST_IN=${EXP}_$(printf "%08d" "$TSTEP")_restart
    ICERST_IN=${OCERST_IN}_ice
    edit_var cn_ocerst_in "'${OCERST_IN}'" ./namelist_cfg
    edit_var cn_icerst_in "'${ICERST_IN}'" ./namelist_ice_cfg

    RST_FILE="none"
    if [ -f $RST_DIR/${OCERST_IN}.nc ] ; then
      RST_FILE=$RST_DIR/${OCERST_IN}.nc
    elif [ -f $RST_DIR/${OCERST_IN}_0000.nc ] ; then
      RST_FILE=$RST_DIR/${OCERST_IN}_0000.nc
    fi
    read -r RST_YR START_MTH START_DAY START_YR <<< "$(date_from_restart $RST_FILE )"
else
    # Start from rest
    edit_var ln_rstart ".false." ./namelist_cfg
    edit_var ln_tsd_init ".true." ./namelist_cfg
    edit_var cn_ocerst_in "'dummy'" ./namelist_cfg
    edit_var cn_icerst_in "'dummy'" ./namelist_ice_cfg
    edit_var nn_it000 "1" ./namelist_cfg
    DATE0=$(get_var nn_date0 ./namelist_cfg)
    START_YR=${DATE0:0:4} ; START_MTH=${DATE0:4:2} ; START_DAY=${DATE0:6:2}
    if [ "$START_MTH" -ne "01" ] ; then
      echo "Warning: not starting at the beginning of the year; itend calculation may not result in an end-of-year finish"
    fi
fi

if [ "$START_YR" -lt 0 ] ; then
    # START_YR not discoverable; assume 365 day year
    YR_LEN=365
else
    YR_LEN=$(greg_length $START_YR )
fi

# Edit nn_itend and (possibly) ln_rst_list and nn_stocklist
# Deal with frequency and length units (convert to hours in all cases)
# TIME_UNITS=0 # 0=years ; 1=days ; 2=hours ; no conversion required for last case.
if [ "$TIME_UNITS" -eq 0 ]; then
    #
    # First calculate the final timestep of the experiment based on LENGTH, ITBEGIN and nn_date0
    # For simplicity, assume always starting at the beginning of the year.
    #
    DATE0=$(get_var nn_date0 ./namelist_cfg)
    EXP_START_YR=${DATE0:0:4} ; EXP_START_MTH=${DATE0:4:2} ; EXP_START_DAY=${DATE0:6:2}
    EXP_LENGTH=$(greg_length $EXP_START_YR )
    LT=$(( LENGTH - 1 ))
    if [ "$LT" -gt 0 ] ; then
      for y in `seq 1 1 $LT`
       do
        EXP_START_YR=$(( EXP_START_YR + 1 ))
        YL=$( greg_length $EXP_START_YR )
        EXP_LENGTH=$(( EXP_LENGTH + YL ))
       done
    fi
    # convert from days to hours
    EXP_LENGTH=$(( EXP_LENGTH * 24 ))
    #
    # Next calculate the end-of-year timesteps in this period of the experiment
    # based on FREQRST, IT000 and ending date of previous job-step
    # Note FREQRST is used to determine the length of each job-step. End of year restarts
    # are always written even when FREQRST > 1
    #
    NEXT_STOCK="$(( IT000 - 1 + YR_LEN * 24 * 60 * 60 / DT ))"
    STOCK_LIST=$NEXT_STOCK
    LT=$(( FREQRST - 1 ))
    if [ "$LT" -gt 0 ] ; then
      for y in `seq 1 1 $LT`
       do
        START_YR=$(( START_YR + 1 ))
        YL=$( greg_length $START_YR )
        NEXT_STOCK="$(( NEXT_STOCK + YL * 24 * 60 * 60 / DT ))"
        # stocklist is limited to length 10
        if [ "$y" -lt 10 ] ; then
          STOCK_LIST="${STOCK_LIST},${NEXT_STOCK}"
        fi
       done
    fi
    ITEND=$NEXT_STOCK
    edit_var nn_stocklist "$STOCK_LIST" ./namelist_cfg
    edit_var ln_rst_list ".true." ./namelist_cfg
elif [ "$TIME_UNITS" -eq 1 ]; then
    ITEND=$(( IT000 - 1 + FREQRST * 24 * 60 * 60 / DT ))
    EXP_LENGTH=$(( LENGTH * 24 ))
    edit_var ln_rst_list ".false." ./namelist_cfg
elif [ "$TIME_UNITS" -eq 2 ]; then
    ITEND=$(( IT000 - 1 + FREQRST * 60 * 60 / DT ))
    EXP_LENGTH=$LENGTH
    edit_var ln_rst_list ".false." ./namelist_cfg
fi

FINAL=$(( ITBEGIN - 1 + EXP_LENGTH * 60 * 60 / DT ))
echo "ITEND = "$ITEND
echo "FINAL = "$FINAL
if [ "$ITEND" -ge "$FINAL" ]; then
    ITEND=$FINAL
fi
edit_var nn_itend "$ITEND" ./namelist_cfg

# Submit next job
if [ "$ITEND" -lt "$FINAL" ] && [ "$RESUBMIT" = true ]; then
    sbatch --dependency=afterok:"$SLURM_JOB_ID" "$SCRIPTNAME"
fi

# Create directories
mkdir -p {RESTARTS,OUTPUT}


# Run checker to trap hanging jobs
./time_step_chk $SLURM_JOB_ID "$SCRIPTNAME" &

srun --overlap --oversubscribe --mem=1GB --ntasks=$SLURM_NNODES --ntasks-per-node=1 ./run_check_mem.sh &
export SLURM_OVERLAP=1

# Run NEMO
echo "Commence srun.."

srun --cpu-bind=verbose --mem-bind=local --overlap --oversubscribe \
--het-group=0 --nodes=12 --ntasks=192 --ntasks-per-core=1 --ntasks-per-node=16 --cpus-per-task=1 /tmp/xios_server.exe \
: --het-group=1 --nodes=34 --ntasks=2924 --ntasks-per-core=1 --ntasks-per-node=86 --cpus-per-task=1  /tmp/nemo \
: --het-group=2 --nodes=1 --ntasks=74 --ntasks-per-core=1 --ntasks-per-node=74 --cpus-per-task=1  /tmp/nemo


# Move log files
mkdir -p LOGS/"$SLURM_JOB_ID"
mv {*.dat,*output*,*.stat*,communication_report.txt} LOGS/"$SLURM_JOB_ID"

# Check time.step and ocean.output
for NAMELIST in *namelist_cfg; do

    # Check that ocean.output exist
    OUTNAME=${NAMELIST//namelist_cfg/ocean.output}
    OUTNAME=LOGS/"$SLURM_JOB_ID"/"$OUTNAME"
    if [ ! -f  "$OUTNAME" ]; then
        echo "E R R O R : $OUTNAME is missing."
        exit 1
    fi

    # Check for errors
    ERROR=$(grep "E R R O R" "$OUTNAME" | wc -l)
    if [ "$ERROR" -gt 0 ]; then
        echo "E R R O R : errors found in $OUTNAME"
        exit 1
    fi

    # Check that time.step exist
    TIMENAME=${NAMELIST//namelist_cfg/time.step}
    if [ ! -f  "$TIMENAME" ]; then
        echo "E R R O R : $TIMENAME is missing."
        exit 1
    fi

    # Check nn_itend vs time.step
    EXPECTED=$(get_var nn_itend "$NAMELIST")
    ACTUAL=$(get_timestep "$TIMENAME")
    if [ "$ACTUAL" -lt "$EXPECTED" ]; then
        echo "E R R O R : $TIMENAME is $ACTUAL (expecting $EXPECTED)"
        exit 1
    fi
done
