import cdsapi
import os
from datetime import date

c = cdsapi.Client()

# Downloads all missing files ERA5 and if latest is true, downloads the latest 3 months in a different directory

#latest = False
latest = True
#first_year = 1940
first_year = 2023
today = date.today()

savedir = '/dssgfs01/scratch/npd/forcing/ERA5/original/'
lastdir = '/dssgfs01/scratch/npd/forcing/ERA5/original_latest/'

# Files needed for NEMO forcing:
variables = ['2m_temperature', '2m_dewpoint_temperature', 
             '10m_u_component_of_wind', '10m_v_component_of_wind',
             'mean_total_precipitation_rate', 'mean_snowfall_rate', 
             'mean_surface_downward_short_wave_radiation_flux', 'mean_surface_downward_long_wave_radiation_flux',
             'mean_sea_level_pressure']

# Files needed for Simon's flux calculations:
#variables = ['2m_temperature', '2m_dewpoint_temperature',
#             '10m_u_component_of_wind', '10m_v_component_of_wind',
#             'mean_total_precipitation_rate', 'mean_sea_level_pressure',
#             'mean_eastward_turbulent_surface_stress', 'mean_northward_turbulent_surface_stress',
#             'mean_surface_latent_heat_flux','mean_surface_sensible_heat_flux',
#             'mean_surface_net_long_wave_radiation_flux', 'mean_surface_net_short_wave_radiation_flux',
#             'sea_ice_cover', 'sea_surface_temperature']

# Files used for wave forcing:
#tier1:
#variables = ['mean_wave_direction', 'mean_wave_period', 'mean_wave_period_based_on_first_moment',
#            'significant_height_of_combined_wind_waves_and_swell',
#            'u_component_stokes_drift', 'v_component_stokes_drift',]
#tier2:
#variables = ['significant_height_of_total_swell', 'significant_height_of_wind_waves',
#            'mean_wave_period_based_on_first_moment_for_swell', 'mean_wave_period_based_on_first_moment_for_wind_waves',
#            'mean_direction_of_total_swell', 'mean_direction_of_wind_waves',
#            'mean_period_of_total_swell','mean_period_of_wind_waves',]
#tier3:
#variables = ['friction_velocity', 'normalized_stress_into_ocean',]

# Compute last complete month:
last_day  = today.day-6  # ERA5 is updated with a 5 day delay, using 6 day offset to make sure and only download completed months
last_year = today.year
last_mon  = today.month

if last_day < 1:
    last_mon = last_mon - 2
else:
    last_mon = last_mon - 1

if last_mon < 1:
    last_year = last_year - 1
    last_mon  = last_mon + 12

save_mon  = last_mon - 3 
save_year = last_year

if save_mon < 1:
    save_year = save_year - 1
    save_mon  = save_mon + 12

if not latest:
  last_mon  = save_mon
  last_year = save_year

for year in range(first_year, last_year + 1):
    if year == last_year:
        lmon = last_mon+1
    else:
        lmon = 13
    for month in range(1, lmon):
        for variable in variables:
            outdir  = (savedir + str(year) + '/' + variable + '/')
            if latest:
                if (((save_year == year) & (save_mon < month)) | (save_year < year)):
                    outdir  = (lastdir + str(year) + '/' + variable + '/')
             
            outfile = (outdir + "{variable}_{year}-{month:02d}.nc".format(variable=variable, year=year, month=month))
            if os.path.isfile(outfile):
                print(outfile + ' exists')
            else:
                if not os.path.exists(outdir):
                    # Create a new directory because it does not exist
                    os.makedirs(outdir)
                print("=========================================================")
                print('Downloading ' +  outfile)
                c.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'variable': variable,
                        'year': str(year),
                        'month': "{month:02d}".format(month=month),
                        'day': [
                            '01', '02', '03',
                            '04', '05', '06',
                            '07', '08', '09',
                            '10', '11', '12',
                            '13', '14', '15',
                            '16', '17', '18',
                            '19', '20', '21',
                            '22', '23', '24',
                            '25', '26', '27',
                            '28', '29', '30',
                            '31',
                        ],
                        'time': [
                            '00:00', '01:00', '02:00',
                            '03:00', '04:00', '05:00',
                            '06:00', '07:00', '08:00',
                            '09:00', '10:00', '11:00',
                            '12:00', '13:00', '14:00',
                            '15:00', '16:00', '17:00',
                            '18:00', '19:00', '20:00',
                            '21:00', '22:00', '23:00',
                        ],
                        'format': 'netcdf',
                    },
                    outfile)
