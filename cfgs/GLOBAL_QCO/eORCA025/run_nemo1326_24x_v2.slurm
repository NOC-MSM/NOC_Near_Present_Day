#!/bin/bash
#SBATCH --job-name=NPD_025
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-socket=32
#SBATCH --nodes=24
##SBATCH --exclude=compute[001-025]

# ========================================================
# PARAMETERS TO SET
# Restart frequency (<0 days, >0 hours)
FREQRST=-1 #1461
# Simulation length (<0 days, >0 hours)
LENGTH=-1   # 16071 == end of 2001
# Parent initial time step (0: infer from time.step)
# IT000 != 0 -> auto-resubmission is switched OFF
IT000=1
# Name of this script (to resubmit)
SCRIPTNAME=run_nemo1326_24x_v2.slurm
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

# Convert days to hours
if [ "$FREQRST" -lt 0 ]; then
    FREQRST=$(( FREQRST * -24 ))
fi
if [ "$LENGTH" -lt 0 ]; then
    LENGTH=$(( LENGTH * -24 ))
fi

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
    OCERST_IN=${EXP}_$(printf "%08d" "$TSTEP")_restart
    ICERST_IN=${OCERST_IN}_ice
    edit_var cn_ocerst_in "'${OCERST_IN}'" ./namelist_cfg
    edit_var cn_icerst_in "'${ICERST_IN}'" ./namelist_ice_cfg

else
    # Start from rest
    edit_var ln_rstart ".false." ./namelist_cfg
    edit_var ln_tsd_init ".true." ./namelist_cfg
    edit_var cn_ocerst_in "'dummy'" ./namelist_cfg
    edit_var cn_icerst_in "'dummy'" ./namelist_ice_cfg
    edit_var nn_it000 "1" ./namelist_cfg
fi

# Edit nn_itend
IT000=$(get_var nn_it000 ./namelist_cfg)
DT=$(get_var rn_Dt ./namelist_cfg)
DT=${DT%.*}
ITEND=$(( IT000 - 1 + FREQRST * 60 * 60 / DT ))
FINAL=$(( LENGTH * 60 * 60 / DT ))
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

# Link nemo
ln -fs ../BLD/bin/nemo.exe nemo

module load NEMO/prg-env

export XIO_HOME=/home/acc/xios-trunk
#
export OMP_NUM_THREADS=1
export I_MPI_SHM=off
unset I_MPI_PMI_LIBRARY 

mpiexec.hydra -print-rank-map -ppn 1 -np 24 $XIO_HOME/bin/xios_server.exe : -np 1326 ./nemo

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

