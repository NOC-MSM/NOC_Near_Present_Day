# **Deep Dives**

!!! abstract "Summary"

    **This page provide additional details on the parameter choices, initial conditions and output options of the Near-Present-Day simulations.**

---

## **Model Configurations :globe_with_meridians:**

The Near-Present-Day integrations have been performed using a traceable hierarchy of three ocean sea-ice model configurations at 1$^{\circ}$, 1/4$^{\circ}$ and 1/12$^{\circ}$ nominal horizontal resolution.

The Near-Present-Day configurations are based on the recent UK Global Ocean and Sea Ice configuration version 9 (GOSI9), which has been documented comprehensively in [Guiavarc’h et al. (2025)](https://doi.org/10.5194/gmd-18-377-2025).

In this section, we describe each of the configurations in detail and provide references for those seeking further technical information.

### **Ocean Model**

All three Near-Present-Day configurations have been performed using the Nucleus for European Modelling of the Ocean ([NEMO](https://www.nemo-ocean.eu)) modelling framework version 4.2.

For a list of the additions made to the version 4.2 release, users are referred to the NEMO [release notes](https://sites.nemo-ocean.io/user-guide/changes.html#changes-since-4-0-7).

#### **Model Grid & Bathymetry**

The ocean model grids and bathymetry used in the Near-Present-Day simulations are the same eORCA global isotropic Mercator grids documented in [Storkey et al (2018)](https://doi.org/10.5194/gmd-11-3187-2018).

The eORCA1, eORCA025 and eORCA12 have a nominal resolution of 1$^{\circ}$, 1/4$^{\circ}$ and 1/12$^{\circ}$ at the Equator. The meridional resolution of the eORCA1 model increases to 1/3$^{\circ}$ at the equator to improve the representation of equatorial ocean dynamics. In the Southern Hemisphere, all grids are extended to 85$^{\circ}$S to include ice shelf cavities (see [Guiavarc’h et al. (2025)](https://doi.org/10.5194/gmd-18-377-2025)).

All three configurations use a non-linear free surface and 75 vertical levels, meaning that grid cell thicknesses throughout the water column vary with time (z$^{*}$ coordinate, [Adcroft and Campin, 2004](https://doi.org/10.1016/j.ocemod.2003.09.003)). Partial step topography is also used to ensure that the depth of each water column equates to the real depth of the ocean by modifying the thickness of the bottom grid cell (e.g., [Barnier et al., 2006](https://doi.org/10.1007/s10236-006-0082-1)).  

#### **Mixing**

* **Momentum Advection:** Vector-invariant form separating horizontal advection into a rotational term (scheme of [Arakawa & Lamb, 1981](https://doi.org/10.1175/1520-0493(1981)109<0018:APEAEC>2.0.CO;2)) and irrotational term (scheme of [Hollingsworth et al., 1983](https://doi.org/10.1002/qj.49710946012)).

* **Tracer Advection:** 4th Order Total Variance Diminishing (TVD) scheme ([Zalesak, 1979](https://doi.org/10.1016/0021-9991(79)90051-2)) in both horizontal and vertical directions.

* **Lateral Diffusion of Momentum:** Performed along geopotential surfaces with the following parameters:

    === "eORCA1"
        * Laplacian viscosity defined by 20,000 m$^{2}$ s$^{-1}$ poleward of 20$^{\circ}$N/S, reducing with the meridional grid size towards the equator.

    === "eORCA025"
        * Bi-Laplacian viscosity with a lateral viscous velocity = 0.0838 m s$^{-1}$

    === "eORCA12"
        * Bi-Laplacian viscosity with a lateral viscous velocity = 0.1895 m s$^{-1}$

* **Lateral Diffusion of Tracers:** Performed on iso-neutral surfaces using Laplacian mixing with the following parameters:

    === "eORCA1"
        * Lateral diffusive velocity = 0.018 m s$^{-1}$.

    === "eORCA025"
        * Lateral diffusive velocity = 0.011 m s$^{-1}$.

    === "eORCA12"
        * Lateral diffusive velocity = 0.027 m s$^{-1}$.

* **Vertical Mixing of Momentum & Tracers:** Performed using an updated version of the [Gaspar et al., 1990](https://doi.org/10.1029/jc095ic09p16179) Turbulent Kinetic Energy (TKE) parameterisation.

* **Adiabatic Eddy Mixing:** Calculated using the parameterisation of [Gent & McWilliams, 1990](https://doi.org/10.1175/1520-0485(1990)020<0150:imiocm>2.0.co;2). 

    === "eORCA1"
        * Uses spatially-dependent coefficients (see [Tréguier et al., 1997](https://doi.org/10.1175/1520-0485(1997)027<0567:POQEIP>2.0.CO;2)).
        * Zonal and meridional eddy-induced velocities (uo_eiv, vo_eiv) are output in addition to the resolved velocities (uo, vo).

    === "eORCA025 & eORCA12"
        * Use spatially- and temporally-dependent version of the parameterisation to account for the unresolved effects of eddies at high latitudes (see [Guiavarc’h et al. (2025)](https://doi.org/10.5194/gmd-18-377-2025)).
        * Zonal and meridional eddy-induced velocities (uo_eiv, vo_eiv) are output in addition to the resolved velocities (uo, vo).

#### **Equation of State**

All three Near-Present-Day simulations use the Thermodynamic Equation Of Seawater 2010 (TEOS-10, [Ioc et al., 2010](https://unesdoc.unesco.org/ark:/48223/pf0000188170)) and hence absolute salinity and conservative temperature variables are output rather than practical salinity and potential temperature as in EOS-80.

#### **Atmospheric Forcing**

Two Near-Present-Day version 1 integrations have been perfomed using the JRA55-do atmospheric reanalysis ([Tsujino et al, 2018](https://doi.org/10.1016/j.ocemod.2018.07.002)) and a climatologically adjusted version of the ERA-5 atmospheric reanalysis ([Hersbach et al., 2020](https://doi.org/10.1002/qj.3803)).

The JRA55-do forced Near-Present-Day integration spans the period 1976-2023, and uses the NCAR bulk formulae ([Large and Yeager, 2009](https://doi.org/10.1007/s00382-008-0441-3)).

The ERA-5 forced Near-Present-Day integration spans the period 1976-present, and uses ECMWF bulk formulae (IFS documentation, cy45). 

There is a well established Surface Air Temperature (SAT) bias at high-latitudes in the ERA-5 atmospheric reanalysis ([Tjernström & Graversen, 2009](https://doi.org/10.1002/qj.380); [Zampieri et al., 2023](https://doi.org/10.1175/MWR-D-22-0130.1)) owing to the poor representation of snow atop of sea ice ([Batrak & Müller, 2019](https://doi.org/10.1038/s41467-019-11975-3)). To account for these large biases, a climatological adjustment is applied to the ERA-5 hourly 2 m temperature field over regions where (ERA-5) sea ice cover > 0%. Climatological offset factors are determined by calculating the difference between the long-term mean (1960-2019) monthly 2 m temperature climatologies of ERA-5 and JRA55-do. To avoid step-like transitions between monthly adjustments, 2 m temperature offset factors are linearly interpolated in time to produce hourly surface forcing fields.

!!! info "Accounting for Current Feedback to the Atmosphere"

    Since the JRA55-do atmospheric reanalysis already takes into account the coupling between ocean currents and surface wind wind stress (termed the Current Feedback), directly forcing an ocean model with reanalysis (relative) winds results in unrealistically weak mesoscale activity and large-scale circulation features (see [Renault et al., 2020](https://doi.org/10.1029/2019MS001715)).

    To overcome this, the parameterisation of [Renault et al., 2017](https://doi.org/10.1038/s41598-017-17939-1) is used to remove the wind and surface stress anomalies induced by the reanalysis surface ocean currents and replace them with those induced by the currents of the Near-Present-Day simulation. This relies on a linear estimate for the current-stress coupling coefficient, $S_{\tau} = \alpha |U_{10_{abs}}| + \beta$, where $\alpha$ = -2.9x10$^{-3}$ N s$^{2}$ m$^{-4}$ and $\beta$ = 0.008 N s m$^{-3}$.

#### **Initial Conditions**

To initialise the Near-Present-Day integrations, conservative temperature and absolute salinity fields from the World Ocean Atlas 2023 ([Reagen et al., 2024](https://www.ncei.noaa.gov/products/world-ocean-atlas)) Climate Normal (30-year average) corresponding to 1971-2000 are used.

Each integration begins with a 3-year spin-up period during which the JRA55-do / adjusted ERA-5 1976 forcing is repeatedly applied while the initial conservative temperature and absolute salinity fields are reset to the final time-step of the previous simulation year.

#### **Sea Surface Salinity Restoration**

The salinity at the sea surface is weakly restored to the World Ocean Atlas 2023 upper 10m average salinity (1991-2020) using a damping magnitude of -33.33 mm day$^{-1}$ in all integrations.

### **Sea Ice Model**

The Near-Present-Day simulations use NEMO's recently introduced dynamic-thermodynamic continuum sea ice model, SI$^{3}$ (Sea Ice modelling Integrated Initiative) (see [Vancoppenolle et al., 2023](https://doi.org/10.5281/zenodo.7534900)).

For a more detailed discussion on the SI$^{3}$ model configuration, users are referred to [Blockley et al. (2023)]() and [Guiavarc’h et al., in review](https://doi.org/10.5194/egusphere-2024-805).

---

## **Creating Initial Conditions :thermometer:**

!!! info "Section Currently Under Development: Come Back Soon!"
    
    **This section will include a description on how to create initial conditions for Near-Present-Day integrations.**

---

## **Editing the Runscript :hammer:**

To run your own Near-Present-Day integrations, you'll likely need to make changes to the ```run_nemo_???.slurm``` runscript in the ```/nemo/cfgs/GLOBAL_QCO/eORCA??/``` directory of your local installation of the NOC_Near_Present_Day repository.

In this section, we discuss the key parameters you will need to modify to begin running your own Near-Present-Day experiments.

Let's start by defining each of the parameters available in an example runscript ```run_nemo_example.slurm```:

```sh title="run_nemo_example.slurm"
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
SCRIPTNAME=run_nemo_example.slurm
# If conducting the repeat and reset T and S spinup set SPIN to 1, else set to 0
SPIN=0
```

* **TIME_UNITS** is used to define the length of the simulation and to specify the resubmission frequency of the runscript. In the above example, ```TIME_UNITS=0``` indicates we are working in years.

* **FREQRST** defines the frequency (in the time units defined above) at which the runscript will be resubmitted as a SLURM batch job. In the above example, ```FREQRST=1``` indicates that every simulation year should be submitted as a separate SLURM batch job (dependent on the successful completion of the previous job/year).

* **IT000** specifies the initial time-step of this job. If ```IT000=0``` then the initial time-step will be determined from the time.step file in the run directory (if no time.step files exists then IT000 is set to 1). When IT000 != 0, runscript resubmission is turned off and its value will be unchanged.

* **ITBEGIN** specifies the starting time-step of the simulation. This is used to update the namelist_cfg file with the final time-step of the job when a restart file is to be written.

* **LENGTH** defines the length of the simulation (in the time units defined above). The simulation length together with the frequency of resubmission (FREQRST) will determine the number of SLURM batch jobs which are created. In the above example, ```LENGTH=3``` and ```FREQRST=1``` means that 3 separate SLURM batch jobs will be required to complete this simulation.

* **SCRIPTNAME** defines the name of the runscript to be resubmitted as a SLURM batch job. In most cases this should be the name of the runscript and should be checked carefully!

* **SPIN** specifies that a spin-up simulation should be perfomed. If ```SPIN=1```, the simulation year will be restarted and the atmospheric forcing repeated, but the initial temperature and salinity fields will be defined from the restart file produced at the final time-step of the previous simulation year. Once a spin-up simulation is completed, define ```SPIN=0``` to continue the simulation with the next year of atmospheric forcing.  

#### A Typical Use Case:

Let's consider a typical example: a user would like to perform a 25-year hindcast simulation (2000-2024) starting with a 5-year spin-up simulation repeating the year 2000. The user would also like each simulation year to be submitted as a separate SLURM batch job.

We can divide this workflow into two runscripts: (1) to perform the spin-up simulation from rest...

```sh title="run_nemo_example_spin-up.slurm"
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
LENGTH=5   
# Name of this script (to resubmit)
SCRIPTNAME=run_nemo_example_spin-up.slurm
# If conducting the repeat and reset T and S spinup set SPIN to 1, else set to 0
SPIN=1
```
...and (2) to perform the 25-year hindcast simulation.

```sh title="run_nemo_example_hindcast.slurm"
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
LENGTH=30   
# Name of this script (to resubmit)
SCRIPTNAME=run_nemo_example_hindcast.slurm
# If conducting the repeat and reset T and S spinup set SPIN to 1, else set to 0
SPIN=0
```

!!! note "Note on Modifying Parameters between Runscripts"
    
    In the example above, we did not modify ```ITBEGIN``` in ```run_nemo_example_hindcast.slurm``` since using ```IT000=0``` means that this simulation will automatically start in 2001 following the final year of the 2000 spin-up simulation.

    Also, note that ```LENGTH=30``` in ```run_nemo_example_hindcast.slurm``` defines the total simulation length in years, including the 5-year spin-up period and 25-year hindcast, rather than the number of additional years to simulate starting from 2001.  

---

## **Storing Output with a Zoom Domain :mag:**

!!! info "Section Currently Under Development: Come Back Soon!"
    
    **This section will include instructions on how to output variables for a sub-domain of the global eORCA grid.**
