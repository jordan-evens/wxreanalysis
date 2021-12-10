#!/usr/bin/env python

import cdsapi
import netCDF4
import datetime
import os
import pandas as pd

LATITUDE_MIN = 41
LATITUDE_MAX = 84
LONGITUDE_MIN = -141
LONGITUDE_MAX = -52
AREA = [LATITUDE_MAX, LONGITUDE_MIN, LATITUDE_MIN, LONGITUDE_MAX]
# NOTE: need to get a client key as per https://cds.climate.copernicus.eu/api-how-to

c = cdsapi.Client()

variable = 'temperature'
year = 2008
hour = 12
csv_file = '{:04d}.csv'.format(year, variable)

def to_datetime(hours):
    return(datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=hours))


MONTHS = list(map(lambda x: '{:02d}'.format(x), range(1, 13)))
DAYS = list(map(lambda x: '{:02d}'.format(x), range(1, 32)))
#
# def get_nc(variable,
#            year,
#            pressure_level=None,
#            months=MONTHS,
#            days=DAYS,
#            hours=[12],
#            level_type='reanalysis-era5-pressure-levels'):
#     nc_file = '{:04d}_{:s}.nc'.format(year,
#                                       '{:04d}'.format(pressure_level) if pressure_level is not None else "sfc")
#     if not os.path.exists(nc_file):
#         args = {
#             'product_type': 'reanalysis',
#             'variable': variable,
#             'year': '{:04d}'.format(year),
#             'month': months,
#             'day': days,
#             'time': list(map(lambda x: '{:02d}:00'.format(x), hours)),
#             'area': AREA,
#             'format': 'netcdf',  # Supported format: grib and netcdf. Default: grib
#         }
#         if pressure_level is not None:
#             args['pressure_level'] = str(pressure_level)
#         print(args)
#         c.retrieve(
#             level_type,
#             args,
#             nc_file)  # Output file. Adapt as you wish.
#     return(nc_file)
area = "/".join(map(str, AREA))
time = list(map(lambda x: '{:02d}:00'.format(x), range(0, 24)))
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'variable': ['temperature', 'relative_humidity'],
        'pressure_level': '2',
        'year': '{:04d}'.format(year),
        'month': MONTHS,
        'day': DAYS,
        'time': time,
        'area': AREA,
        'format': 'netcdf',
    },
    '{year:04d}_temp_rh.nc'.format(year=year))
# get_nc('temperature', 2008, 2)
# get_nc('relative_humidity', 2008, 2)
c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'variable': ['u_component_of_wind', 'v_component_of_wind'],
        'pressure_level': '10',
        'year': '{:04d}'.format(year),
        'month': MONTHS,
        'day': DAYS,
        'time': time,
        'area': AREA,
        'format': 'netcdf',
    },
    '{year:04d}_wind.nc'.format(year=year))
# get_nc('u_component_of_wind', 2008, 10)
# get_nc('v_component_of_wind', 2008, 10)
c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': 'total_precipitation',
        'year': '{:04d}'.format(year),
        'month': MONTHS,
        'day': DAYS,
        'time': time,
        'area': AREA,
        'format': 'netcdf',
    },
    '{year:04d}_prec.nc'.format(year=year))
#get_nc('total_precipitation', 2008, hours=range(0, 24), level_type='reanalysis-era5-single-levels')
#
# data = netCDF4.Dataset(nc_file)
# start_date = to_datetime(int(data.variables['time'][[0]].data[0]))
# key = list(data.variables.keys())[3]
# df = pd.DataFrame(data=None, columns=["date", "latitude", "longitude", key])
# df.to_csv(csv_file, index=False)
# #result = []
# for i_time, time in enumerate(data.variables['time']):
#     date = to_datetime(int(time.data))
#     print(date)
#     by_latitude = []
#     for i_lat, latitude in enumerate(data.variables['latitude']):
#         longitudes = data.variables['longitude'][:].data
#         values = data.variables[key][i_time, i_lat].data
#         df = pd.DataFrame({'longitude': longitudes, variable: values})
#         df['latitude'] = latitude
#         by_latitude.append(df)
#     df = pd.concat(by_latitude)
#     df['time'] = date
#     #result.append(df)
#     df = df[['date', 'latitude', 'longitude', variable]]
#     df.to_csv(csv_file, header=False, index=False, mode='a')
#
