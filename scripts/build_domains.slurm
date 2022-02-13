#!/bin/bash
#SBATCH --job-name=DOMAINcfg
#SBATCH --time=00:20:00
#SBATCH --account=n01-CLASS
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --ntasks=128

# Set environment
module swap PrgEnv-cray/8.0.0 PrgEnv-gnu/8.1.0
module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.0.7
module load cray-netcdf-hdf5parallel/4.7.4.7
module load libfabric
module list
export OMP_NUM_THREADS=1

# Link exe
ln -fs ../make_domain_cfg.exe .
ln -fs ../../REBUILD_NEMO/rebuild_nemo{.exe,} .

# Create domains
srun --ntasks="$SLURM_NTASKS" ./make_domain_cfg.exe || exit

# Rebuild
for FILE in *namelist_cfg; do
    STRING=$(basename -- "$FILE")
    PATTERN="domain_cfg"
    PREFIX="${STRING/namelist_cfg/$PATTERN}"
    if ./rebuild_nemo -n nam_rebuild -p "$SLURM_NTASKS" -d 1 -x 64 -y 64 -z 1 -t 1 "$PREFIX" "$SLURM_NTASKS"; then
        rm "${PREFIX}"_[0-9][0-9][0-9][0-9].nc
    else
        rm "$PREFIX".nc
    fi
done

# Bye-Bye
sacct -j "$SLURM_JOB_ID" --units=G --format=JobID,JobName,NNodes,Elapsed,AveRSS,MaxRSS,MaxRSSNode,MaxRSSTask
