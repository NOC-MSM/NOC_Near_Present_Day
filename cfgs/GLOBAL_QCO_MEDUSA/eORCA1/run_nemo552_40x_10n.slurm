#!/bin/bash
#SBATCH --job-name=NPD_1
#SBATCH --partition=test
#SBATCH --time=00:20:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=10
##SBATCH --exclude=compute[001-059]
#SBATCH --exclude=compute[001]

# ========================================================
# PARAMETERS TO SET
# time units used here for restart frequency and simulaion length
  TIME_UNITS=0 # 0=years ; 1=days ; 2=hours
# Restart/resubmission frequency (in TIME_UNITS)
  FREQRST=1
# job-step initial time step (0: infer from time.step)
# IT000 != 0 -> auto-resubmission is switched OFF
  IT000=0
#
# Simulation original starting time step (unchanged for LENGTHxTIME_UNITS)
  ITBEGIN=1
# Simulation length (in TIME_UNITS) 
  LENGTH=10
# Name of this script (to resubmit)
  SCRIPTNAME=run_nemo552_40x_10n.slurm
# If conducting the repeat and reset T and S spinup set SPIN to 1, else set to 0
  SPIN=0 
# Move OUTPUTS to a store-place :                            
  ISTORE=.false.                                               
export STORE_OUT="/general/n01/n01/julpal/OUTPUTS/NPD/eORCA025_river"  
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
module purge
module load NEMO/prg-env
module load libfabric/1.18.0-GCCcore-12.3.0

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

greg_month_length() {
    # Returns month length (in days) given month and year (Gregorian calendar)
    MN=$1
    YR=$2
    if [ $MN -eq 1 ] ; then  
      echo 31
    elif [ $MN -eq 2 ] ; then  
      if [ $(( $YR % 400 )) -eq 0 ] ; then
        echo 29
      elif [ $(( $YR % 100 )) -eq 0 ] ; then
        echo 28
      elif [ $(( $YR % 4 )) -eq 0 ] ; then
        echo 29
      else
        echo 28
      fi
    elif [ $MN -eq 3 ] ; then  
      echo 31
    elif [ $MN -eq 4 ] ; then  
      echo 30
    elif [ $MN -eq 5 ] ; then  
      echo 31
    elif [ $MN -eq 6 ] ; then  
      echo 30
    elif [ $MN -eq 7 ] ; then  
      echo 31
    elif [ $MN -eq 8 ] ; then  
      echo 31
    elif [ $MN -eq 9 ] ; then  
      echo 30
    elif [ $MN -eq 10 ] ; then  
      echo 31
    elif [ $MN -eq 11 ] ; then  
      echo 30
    elif [ $MN -eq 12 ] ; then  
      echo 31
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

## END useful Functions
# ========================================================
# Edit namelists

if [ "$IT000" -gt 1 ]; then
    # Use restarts
    edit_var ln_rstart ".true." ./namelist_cfg
    edit_var ln_tsd_init ".false." ./namelist_cfg
    edit_var ln_rsttr ".true." ./namelist_top_cfg
    edit_var ln_trcdta ".false." ./namelist_top_cfg

    # Set start
    TSTEP=$((IT000 - 1))
    IT000=$((TSTEP + 1))
    edit_var "nn_it000" "$IT000" ./namelist_cfg

    # Set restart prefixes
    EXP=$(get_var cn_exp ./namelist_cfg)
    RST_DIR=$(get_var cn_ocerst_indir ./namelist_cfg)
    OCERST_IN=${EXP}_$(printf "%08d" "$TSTEP")_restart
    ICERST_IN=${OCERST_IN}_ice
    TCERST_IN=${OCERST_IN}_trc
    edit_var cn_ocerst_in "'${OCERST_IN}'" ./namelist_cfg
    edit_var cn_icerst_in "'${ICERST_IN}'" ./namelist_ice_cfg
    edit_var cn_trcrst_in "'${TCERST_IN}'" ./namelist_top_cfg

    RST_FILE="none"
    if [ -f $RST_DIR/${OCERST_IN}.nc ] ; then
      RST_FILE=$RST_DIR/${OCERST_IN}.nc
    elif [ -f $RST_DIR/${OCERST_IN}_0000.nc ] ; then
      RST_FILE=$RST_DIR/${OCERST_IN}_0000.nc
    fi
    read -r RST_YR START_MTH START_DAY START_YR <<< "$(date_from_restart $RST_FILE )"

    if [ $SPIN = 1 ] ; then
      START_YR=$RST_YR
      edit_var ln_reset_ts ".true." ./namelist_cfg
      edit_var nn_rstctl "0" ./namelist_cfg
    elif [ $SPIN = 0 ] ; then
      edit_var ln_reset_ts ".false." ./namelist_cfg
      edit_var nn_rstctl "2" ./namelist_cfg
    else
      echo "Error: it is unclear whether this is a spinup with T and S reset, or a normal simulation"
      exit 1
    fi

    echo "Start year is: " $START_YR

else
    # Start from rest
    edit_var ln_rstart ".false." ./namelist_cfg
    edit_var ln_tsd_init ".true." ./namelist_cfg
    edit_var ln_rsttr ".false." ./namelist_top_cfg
    edit_var ln_trcdta ".true." ./namelist_top_cfg
    edit_var cn_ocerst_in "'dummy'" ./namelist_cfg
    edit_var cn_icerst_in "'dummy'" ./namelist_ice_cfg
    edit_var cn_trcrst_in "'dummy'" ./namelist_top_cfg
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
      if [ $SPIN = 1 ] ; then
        EXP_LENGTH=$(( EXP_LENGTH * LENGTH ))
      else
        for y in `seq 1 1 $LT`
        do
          EXP_START_YR=$(( EXP_START_YR + 1 ))
          YL=$( greg_length $EXP_START_YR )
          EXP_LENGTH=$(( EXP_LENGTH + YL ))
        done
      fi
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
    if [ "$ISTORE" ]; then
      sbatch --dependency=afterok:"$SLURM_JOB_ID" "e025tordf.slurm" ${STORE_OUT}
    fi
fi

# Create directories
mkdir -p {RESTARTS,OUTPUT}

# Link nemo
ln -fs ../BLD/bin/nemo.exe nemo

export XIO_HOME=/dssgfs01/working/atb299/XIOS3
#
export OMP_NUM_THREADS=1
export I_MPI_SHM=off
unset I_MPI_PMI_LIBRARY 

# Run checker to trap hanging jobs
./time_step_chk $SLURM_JOB_ID "$SCRIPTNAME" &

mpiexec.hydra -print-rank-map -ppn 1 -np 57 $XIO_HOME/bin/xios_server.exe : -np 552 ./nemo

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

