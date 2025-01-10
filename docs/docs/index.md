# Getting Started

**Welcome to the documentation for the Near-Present Day simulations developed by the National Oceanography Centre :wave:**

## Introduction :ocean:
As part of the Atlantic Climate and Environment Strategic Science ([**AtlantiS**](https://noc.ac.uk/projects/atlantis)) project, the National Oceanography Centre is developing a full suite of global ocean model configurations to perform multi-decadal Near-Present-Day simulations.

Our aim is that the Near-Present-Day simulations will be kept up to date with a 1-3 month lag.

---

## Configurations :globe_with_meridians:

The Near-Present-Day simulations consists of a hierarchy of three ocean sea-ice configurations of NEMO v4.2 at 1$^{\circ}$, 1/4$^{\circ}$ and 1/12$^{\circ}$ nominal horizontal resolution.

The key features of each configuration are summarised below:

=== "eORCA1"
    * 1$^{\circ}$ nominal horizontal resolution (j=331, i=360).
    * 75 vertical z$^{*}$ levels.
    * Eddy induced velocities determined using the [Gent and McWilliams (1990)](https://doi.org/10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2) diffusion scheme.
    * Coupled to [SI$^{3}$](https://doi.org/10.5281/zenodo.7534900) sea ice engine.
    * Initialised from [World Ocean Atlas 2023](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/) (1971-2001) climatology.
    * Forced with JRA55-do (v1; 1976-2023) and climatologically corrected ERA-5 (v2; 1976-present) atmospheric forcing.

=== "eORCA025"
    * 1/4$^{\circ}$ nominal horizontal resolution (j=1206, i=1440).
    * 75 vertical z$^{*}$ levels.
    * Eddy induced velocities determined using the grid-scale dependent [Gent and McWilliams (1990)](https://doi.org/10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2) diffusion scheme.
    * Coupled to [SI$^{3}$](https://doi.org/10.5281/zenodo.7534900) sea ice engine.
    * Initialised from [World Ocean Atlas 2023](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/) (1971-2001) climatology.
    * Forced with JRA55-do (v1; 1976-2023) and climatologically corrected ERA-5 (v2; 1976-present) atmospheric forcing.

=== "eORCA12"
    * 1/12$^{\circ}$ nominal horizontal resolution (j=3605, i=4320).
    * 75 vertical z$^{*}$ levels.
    * Coupled to [SI$^{3}$](https://doi.org/10.5281/zenodo.7534900) sea ice engine.
    * Initialised from [World Ocean Atlas 2023](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/) (1971-2001) climatology.
    * Forced with JRA55-do (v1; 1976-2023) and climatologically corrected ERA-5 (v2; 1976-present) atmospheric forcing.

For more details on each model configuration see [Deep Dives: Model Configurations].

[Deep Dives: Model Configurations]: deep_dives.md#model-configurations

---

## Quick Start :rocket:

### Installation

To get started, check out and set up an instance of the NPD GitHub [repository](https://github.com/NOC-MSM/NOC_Near_Present_Day):

```sh
git clone git@github.com:NOC-MSM/NOC_Near_Present_Day.git
```

??? tip "Helpful Tip..."

    * **It is not advised to checkout the respository in your home directory.**

Next, run the setup script to download [NEMO](https://www.nemo-ocean.eu) & compile the tools and configurations:

=== "Anemone"
    ```sh
    cd NOC_Near_Present_Day

    ./setup
    ```

=== "Archer2"
    ```sh
    cd NOC_Near_Present_Day

    ./setup -s Archer2
    ```

By default, the script will setup for the Anemone HPC, which is ideally suited for development tasks or for the fast turnaround of smaller NPD configurations (e.g., eORCA1 or eORCA025).

### Running An Experiment

The global eORCA1 and eORCA025 configurations are ready to run. Here, we provide a brief overview on how to setup a first experiment with the default atmospheric forcing JRA55-do on the Anemone HPC:

=== "eORCA1"
    All that is required to run the eORCA1 Near-Present-Day configuration starting from 1976 is:
    ```sh
    cd nemo/cfgs/GLOBAL_QCO/eORCA1

    sbatch run_nemo552_40x_10n.slurm
    ```

    There are a few important variables to set in the runscript.
    
    The default configuration will generate a 10-year simulation without spin-up which is divided into 1-year jobs:
    ```sh
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
    ```

    To run with ERA-5 atmospheric forcing, use the following command to modify the namelist file:
    ```sh
    cp namelist_cfg_ERA5 namelist_cfg
    ```

=== "eORCA025"
    All that is required to run the eORCA025 Near-Present-Day configuration starting from 1976 is:
    ```sh
    cd nemo/cfgs/GLOBAL_QCO/eORCA025

    sbatch run_nemo1326_96x_26n.slurm
    ```

    There are a few important variables to set in the runscript.
    
    The default configuration will generate a 3-year simulation without spin-up which is divided into 1-year jobs:
    ```sh
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
    LENGTH=3   
    # Name of this script (to resubmit)
    SCRIPTNAME=run_nemo1326_96x_26n.slurm
    # If conducting the repeat and reset T and S spinup set SPIN to 1, else set to 0
    SPIN=0
    ```

    To run with ERA-5 atmospheric forcing, use the following command to modify the namelist file:
    ```sh
    cp namelist_cfg_ERA5 namelist_cfg
    ```



