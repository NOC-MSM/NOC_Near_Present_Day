# NOC Near-Present Day: OSNAP Diagnostics

## Last Edited:
**Ollie Tooth (oliver.tooth@noc.ac.uk) |  06-08-2025**

## Description:
Instructions to calculate offline diagnostics along the Overturning in the Subpolar North Atlantic Program (**OSNAP**) hydrographic section using monthly mean outputs from the National Oceanography Centre Near-Present Day simualation.

### Step 1. Preparing Python Virtual Environment

- To extract the OSNAP section from NEMO ocean model outputs, we will use the `nemo_cookbook` Python package, which should be cloned from GitHub as follows:

```bash
git clone git@github.com:NOC-MSM/nemo_cookbook.git
```

- Create a new Python virtual environment using miniforge:

```bash
conda create -n env_npd_diags pyhon=3.13
```

- Locally install the nemo_cookbook Python package from the cloned repository:

```bash 
cd nemo_cookbook

pip install -e .
```

### Step 2. Calculating OSNAP Offline Diagnostics

- **Anemone Users:** Run the eORCA1 & eORCA025 bash script `NOC_Near_Present_Day/scripts/diagnostics/osnap/eORCA1/run_extract_osnap_eORCA1.sh` to submit a SLURM job to the compute queue to perform the OSNAP extraction in parallel using Dask and output to a netCDF file.

```bash
cd /NOC_Near_Present_Day/scripts/diagnostics/osnap/eORCA1

./run_extract_osnap_eORCA1.sh
```

- **Archer2 Users:** Run the eORCA12 bash script `NOC_Near_Present_Day/scripts/diagnostics/osnap/eORCA1/run_extract_osnap_eORCA1.sh` to submit a SLURM job to the Data Analysis queue perform the OSNAP extraction in parallel using Dask and output to a netCDF file.

```bash
cd /NOC_Near_Present_Day/scripts/diagnostics/osnap/eORCA12

./run_extract_osnap_eORCA12.sh
```