# NOC Initial Conditions for NEMO eORCA simulations

Generating initial conservative temperature and absolute salinity fields and sea surface salinity restoration conditions for global configurations at eORCA1, eORCA025 and eORCA12 using World Ocean Atlas climate norms.

## Quick start on {Archer2|Anemone}

The following commands will check out and set up an instance of NPD. It is not advised to do this in your home directory.

```shell
git clone git@github.com:NOC-MSM/NOC_Near_Present_Day.git
cd NOC_Near_Present_Day
./setup {-s Archer2}
```

## Downloading World Ocean Atlas 23 climate norms

First, we can use the following commands to run a bash script which will download the upper 1500m monthly mean climatologies for the climate norms (1971-2000, 1981-2010, 1991-2020) and the winter-average full-depth climatologies from the [NCEI](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/).

```shell
cd tools/INITIAL_CONDITIONS/

./download_woa23_climatologies.sh
```

## Setting Up a Virtual Environment

Next, use the following commands to set-up the conda virtual environment needed to run the initial condition and SSS restoration Python scripts.

```shell
cd tools/INITIAL_CONDITIONS/

conda env create -f env_initial_conditions.yml

conda activate env_initial_conditions
```

## Creating Initial Conditions

```create_initial_conditions.py``` produces monthly mean initial conservative temperature (t_con) and absolute salinity (s_abs) fields by combining a chosen WOA23 climate norm in the upper 1500m with the winter-average properties for the full ocean depth. The output is saved to a single netCDF file in the working directory.

Once you have modified the user inputs included at the top of the script, this can be run using your virtual environment as follows:

```shell
python create_initial_conditions.py
```

## Creating SSS Restoring Conditions

```create_sss_restoring_climatology.py``` produces monthly mean upper 10m absolute salinity (s_abs) fields using a chosen upper 1500m WOA23 climate norm. The output is saved to a single netCDF file in the working directory.

Once you have modified the user inputs included at the top of the script, this can be run using your virtual environment as follows:

```shell
python create_sss_restoring_climatologies.py
```

**Note**: In both scripts, users have the option of reading the 3-dimensional indexes used to flood fill grid cells on land with the values of a variable at the nearest ocean grid cell from a netCDF file or generating them and writing to a netCDF file for future use. This process is not dependent on your target model's horizontal resolution (i.e., the same flood fill indexes can be used with eORCA1 and eORCA12 as they apply only the WOA grid and are used prior to regridding.)