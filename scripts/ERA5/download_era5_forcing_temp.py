import cdsapi
import os
from datetime import date

c = cdsapi.Client()

first_year = 1940
today = date.today()
last_year = today.year
last_mon  = today.month-4
if last_mon < 1:
    last_year = last_year - 1
    last_mon  = last_mon + 12

variables = ['2m_temperature', '2m_dewpoint_temperature', 
             '10m_u_component_of_wind', '10m_v_component_of_wind',
             'mean_total_precipitation_rate', 'mean_snowfall_rate', 
             'mean_surface_downward_short_wave_radiation_flux', 'mean_surface_downward_long_wave_radiation_flux',
             'mean_sea_level_pressure']

for year in range(first_year, last_year + 1):
    if year == last_year:
        lmon = last_mon+1
    else:
        lmon = 13 
    for month in range(6, lmon):
        for variable in variables:
            outdir  = ('/dssgfs01/scratch/npd/forcing/ERA5/original/' + str(year) + '/' + variable + '/')
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
