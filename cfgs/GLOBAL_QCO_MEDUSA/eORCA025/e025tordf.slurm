#!/bin/bash

# Slurm job options (job-name, job time)
#SBATCH --job-name=archit
#SBATCH --time=2:59:0
#SBATCH --ntasks=1

#SBATCH --account=n01-CLASS
#SBATCH --partition=serial
#SBATCH --qos=serial

# Define memory required for this jobs. By default, you would
# get just under 2 GB, but you can ask for up to 125 GB.
#SBATCH --mem=8G

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
#
srcdir=${PWD}/OUTPUT/
dstroot=$1  ## target storing directory
ecount=0
cd $srcdir

## fresh kill list
rm to_be_deleted.txt
#
echo ptrc files to proceed : 
ls *ptrc*nc
#
for year in `seq 1976 1 2023`
do
 dstdir=${dstroot}/${year}
#
 if [ ! -d $dstdir ] ; then mkdir -p $dstdir ; fi
#
 for t in 1y 1m 5d 1d
 do
  for g in scalar icemod grid_U grid_V grid_W grid_T diad_T ptrc_T
  do
   files=`ls -1 *${t}_${g}_${year}* 2>/dev/null`
   #files=`ls \*${t}_${g}_${year}\* `
   #echo ${t}_${g}_${year}
   if [ ${#files} -lt 1 ] ; then
    echo 'Nothing to transfer for '$g' at '${t}' interval in '${year}
   else
    scount=0
    numfiles=`wc -w <<< ${files}`
    echo 'Transferring '${numfiles}' for '$g' at '${t}' interval:'
    for f in $files
    do
     if [ ! -f $dstdir/$f ] ; then
      echo rsync -vtpo $f $dstdir
      rsync -vtpo  $f $dstdir
      #cp -p $f $dstdir
      if [ $? -eq 0 ] ; then scount=$(( $scount + 1 )) ; fi
     else
      ## file already there check if files are the same size : 
      heres=`wc -c $f | cut -d' ' -f1`
      theres=`wc -c $dstdir/$f | cut -d' ' -f1`
      if [ $heres -eq $theres ]
        then
        #echo $f -- copy size match -- can safely erase local copy
        echo rm -f $f >> to_be_deleted.txt
        # rm -f $f
      else 
        ## file broken -- try to complete it 
        echo fixing broken copy of $f with
        echo rsync -vtpo --append $f $dstdir
        rsync -vtpo --append $f $dstdir
      fi 
     fi
    done
    if [ $scount -ne $numfiles ] ; then echo 'Errors have occurred' ; ecount=$(( $ecount + 1 )) ; fi
   fi
  done
 done
done
#
if [ $ecount -gt 0 ] ; then 
 echo 'ERRORS HAVE OCCURRED. Check log for Errors'
else
 echo 'Transfers successful'
fi


# Find and delete all files of zero size in the current directory and subdirectories
cd ${dstroot}
echo 'Erasing empty copies --'
find . -type f -size 0 -print -delete




