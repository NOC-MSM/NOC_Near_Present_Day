# NOC Near-Present-Day simulation

## Quick start on {Archer2|Anemone} (UK National Supercomputing Service)
```shell
git clone git@github.com:NOC-MSM/NOC_Near_Present_Day.git -b NEMO_v4.2
cd NOC_Near_Present_Day
./setup -s {Archer2|Anemone}
```
The setup script downloads nemo, compiles tools and configurations. Note the "--gnu" option may be necessary, depending on system and compiler choice. 

To run NEMO:
```shell
cd nemo/cfgs/GLOBAL_QCO/eORCA12
```
and create a runscript.

Example `mkslurm_NPD` settings for production runs on Archer2:
```shell
../../../scripts/python/mkslurm_NPD -S 48 -s 16 -m 1 -C 5504 -g 0 -a n01-CLASS -j eORCA12 -t 1-00:00:00 --gnu > run_nemo5504_48X.slurm

../../../scripts/python/mkslurm_NPD -S 48 -s 16 -m 1 -C 8238 -g 0 -a n01-CLASS -j eORCA12 -t 0-00:10:00 --gnu > run_nemo8238_48X.slurm

../../../scripts/python/mkslurm_NPD -S 48 -s 16 -m 1 -C 11168 -g 0 -a n01-CLASS -j eORCA12 -t 0-00:10:00 --gnu > run_nemo11168_48X.slurm
```

Example `mkslurm_NPD` settings for production runs on Anemone:
```shell
../../../scripts/python/mkslurm_NPD -S 24 -m 1 -C 3496 -j eORCA12 -t 00:20:00 > run_nemo3496_24X.slurm
```


There are a few variables to set in `run_nemo.slurm`. For example, the following variables will generate a 2-hour simulation split in 1-hour jobs.
```bash
# ========================================================
# PARAMETERS TO SET
# Restart frequency (<0 days, >0 hours)
FREQRST=1
# Simulation length (<0 days, >0 hours)
LENGTH=2
# Parent initial time step (0: infer from time.step)
# PARENT_IT000 != 0 -> auto-resubmission is switched OFF
PARENT_IT000=0
# Name of this script (to resubmit)
SCRIPTNAME=run_nemo-short.slurm
# =======================================================
```
Finally:
```shell
sbatch run_nemo.slurm
```




## Setup
### Global eORCA12
Resolution:
- Horizontal: 1/12°
- Vertical: 75 levels
