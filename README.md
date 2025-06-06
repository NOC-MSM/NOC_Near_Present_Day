<br />
<p align="center">
    <img src="./docs/docs/assets/icons/noc_logo.png" alt="Logo" width="240" height="80">
  </a>


  <h1 align="center">Near-Present-Day Simulations</h1>

  <p align="center">
    Global ocean model configurations (eORCA1, eORCA025 & eORCA12) to perform multi-decadal Near-Present-Day simulations.
    </a>
    <br />
    <br />
    ·
    <a href="https://noc-msm.github.io/NOC_Near_Present_Day/"><strong>Read the Docs</strong></a>
    ·
    <a href="https://github.com/NOC-MSM/NOC_Near_Present_Day/issues"><strong>Report an Issue</strong></a>
    ·
    <br />
    <br />
    <a href="https://doi.org/10.5281/zenodo.15310353"><img src="https://zenodo.org/badge/458888476.svg" alt="DOI"></a>
  </p>
</p>


## Quick start on {Archer2|Anemone}

The following commands will check out and set up an instance of NPD. It is not advised to do this in your home directory.

```shell
git clone git@github.com:NOC-MSM/NOC_Near_Present_Day.git

cd NOC_Near_Present_Day

./setup {-s Archer2}
```
The setup script downloads nemo, compiles tools and configurations. Setup defaults to Anemone, which is ideally suited for fast development/turnaround of smaller configurations (e.g. eORCA025). 

### Global eORCA025

The global eORCA025 configuration is ready to run. All that is required is:

```shell
cd nemo/cfgs/GLOBAL_QCO/eORCA025

sbatch run_nemo1326_24x_v2.slurm
```

There are a few variables to set in the runscript. For example, the following variables will generate a 2-hour simulation split in 1-hour jobs.
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
SCRIPTNAME=run_nemo.slurm
# =======================================================
```

The default forcing is currently JRA55-do. To run with ERA5 instead:

```
cp namelist_cfg_ERA5 namelist_cfg
```

### Global eORCA12

The target system for Global eORCA12 should be Archer2. It is possible to run on Anemone, but not recommended.

To run NEMO:
```shell
cd nemo/cfgs/GLOBAL_QCO/eORCA12
```
and create a runscript.

Example `mkslurm_NPD` settings for eORCA12 production runs on Archer2:
```shell
../../../scripts/python/mkslurm_NPD -S 48 -s 16 -m 1 -C 5504 -g 0 -a n01-CLASS -j eORCA12 -t 1-00:00:00 --gnu > run_nemo5504_48X.slurm

../../../scripts/python/mkslurm_NPD -S 48 -s 16 -m 1 -C 8238 -g 0 -a n01-CLASS -j eORCA12 -t 0-00:10:00 --gnu > run_nemo8238_48X.slurm

../../../scripts/python/mkslurm_NPD -S 48 -s 16 -m 1 -C 11168 -g 0 -a n01-CLASS -j eORCA12 -t 0-00:10:00 --gnu > run_nemo11168_48X.slurm
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
sbatch run_nemo11168_48X.slurm
```

## Using Cylc Workflow [Under development]

NPD can be run using Cylc, which provides the ability to graph the workflow, automating the resubmission of job increments, post-processing, data transfer, etc. It is currently being set up and tested on Anemone. The NEMO_cylc workflow is accessible [here](https://github.com/NOC-OI/NEMO_cylc).

## Current Configurations
### Global eORCA1
Resolution:
- Horizontal: 1°
- Vertical: 75 levels

### Global eORCA025
Resolution:
- Horizontal: 1/4°
- Vertical: 75 levels

### Global eORCA12
Resolution:
- Horizontal: 1/2°
- Vertical: 75 levels

For more information see the [Near-Present-Day Documentation](https://noc-msm.github.io/NOC_Near_Present_Day/).
