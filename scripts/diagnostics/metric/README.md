# NOC Near-Present Day: METRIC Diagnostics

## Last Edited:
**Ollie Tooth (oliver.tooth@noc.ac.uk) |  06-08-2025**

## Description
Instructions to calculate offline diagnostics along the RAPID-MOCHA 26.5N, MOVE 16N and SAMBA 34.5S hydrographic sections using monthly mean outputs from the National Oceanography Centre Near-Present Day simualation.

### Step 1. Preparing Python Virtual Environment

- To extract sections from NEMO ocean model outputs, we will use a forked version of the `metric` Python package which uses xarray for opening multi-file datasets. The library should be cloned from GitHub as follows:

```bash
cd /NOC_Near_Present_Day/scripts/diagnostics/metric/
git clone git@github.com:oj-tooth/metric.git
```

- Create a new Python virtual environment using miniforge:

```bash
conda create -n env_npd_diags pyhon=3.13
```

- Locally install the metric Python package from the cloned repository:

```bash
cd metric

pip install -e .
```

### Step 2. Calculating OSNAP Offline Diagnostics

- **Anemone Users:** First, modify Line 39 in each of the eORCA1 & eORCA025 bash scripts in `NOC_Near_Present_Day/scripts/diagnostics/metric/eORCA1/run_metric_*.slurm` and `NOC_Near_Present_Day/scripts/diagnostics/metric/eORCA025/run_metric_*.slurm`. Then, submit SLURM jobs to the compute queue to produce observational-equivalent diagnostics for each section in parallel using Dask and output to a netCDF file.

```bash
cd /NOC_Near_Present_Day/scripts/diagnostics/metric/eORCA1

./run_metric_RAPID_eORCA1_ERA5_v1.slurm
```

- **Archer2 Users:** Modify Line 39 of the eORCA12 bash scripts `NOC_Near_Present_Day/scripts/diagnostics/metric/eORCA1/run_metric_*.slurm` and submit SLURM jobs to the Data Analysis queue to produce observational-equivalent diagnostics for each section in parallel using Dask and output to a netCDF file.

```bash
cd /NOC_Near_Present_Day/scripts/diagnostics/osnap/eORCA12

./run_metric_RAPID_eORCA12_ERA5_v1.slurm
```

## Note:

- See the `config/` directory for Jupyter Notebooks used to define the RAPID-MOCHA 26.5N, MOVE 16N and SAMBA 34.5N array indexes on each NEMO model grid.