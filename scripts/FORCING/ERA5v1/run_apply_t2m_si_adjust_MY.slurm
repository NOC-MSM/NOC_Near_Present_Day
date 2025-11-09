#!/bin/bash
#SBATCH --job-name=t2m_si_adjust_MY
#SBATCH --time=01:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --exclusive

# ------------------------------------------------------------
# run_apply_t2m_si_adjust_MY.slurm
# Description: Apply climatological adjustment to ERA-5 2m air
# temperature over sea ice for multiple years simultaneously.
#
# Created By: Adam Blaker (atb299@noc.ac.uk)
# Last Edited By: Ollie Tooth (oliver.tooth@noc.ac.uk)
# ------------------------------------------------------------
# -- Load Dependencies -- #
module load NEMO/prg-env
export I_MPI_SHM=icx

source /home/otooth/miniconda3/etc/profile.d/conda.sh
conda activate /dssgfs01/working/otooth/conda_envs/env_npd_era5

# -- Process -P N Years Simultaneously -- #
# Define years to process:
yrs=($(seq 2023 2025))

echo "================================================="
echo "In Progress: Applying Climatological Adjustment to ${yrs} ..."
xargs -n 1 -P 2 <<<${yrs[*]} sh -c 'python3 apply_ERA5_t2m_si_adjustment.py --year ${1}' bash
echo "Completed: Applied Climatological Adjustment to ${yrs}."
