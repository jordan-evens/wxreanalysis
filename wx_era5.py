#!/usr/bin/env python

import ssl

## So HTTPS transfers work properly
ssl._create_default_https_context = ssl._create_unverified_context

import requests
from urllib3.exceptions import InsecureRequestWarning

# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)

import cdsapi
import netCDF4
import datetime
import os
import pandas as pd
import numpy as np
import math

LATITUDE_MIN = 41
LATITUDE_MAX = 84
LONGITUDE_MIN = -141
LONGITUDE_MAX = -52
AREA = [LATITUDE_MAX, LONGITUDE_MIN, LATITUDE_MIN, LONGITUDE_MAX]
# NOTE: need to get a client key as per https://cds.climate.copernicus.eu/api-how-to

def ensure_dir(dir):
    """!
    Check if directory exists and make it if not
    @param dir Directory to ensure existence of
    @return None
    """
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir

DIR = ensure_dir('era5')

def calc_rh(Td, T):
    """!
    Calculate RH based on dewpoint temp and temp
    @param Td Dewpoint temperature (Celcius)
    @param T Temperature (Celcius)
    @return Relative Humidity (%)
    """
    m = 7.591386
    Tn = 240.7263
    return 100 * pow(10, (m * ((Td / (Td + Tn)) - (T / (T + Tn)))))


def calc_ws(u, v):
    """!
    Calculate wind speed in km/h from U and V wind given in m/s
    @param u Wind vector in U direction (m/s)
    @param v Wind vector in V direction (m/s)
    @return Calculated wind speed (km/h)
    """
    # NOTE: convert m/s to km/h
    return 3.6 * math.sqrt(u * u + v * v)


def calc_wd(u, v):
    """!
    Calculate wind direction from U and V wind given in m/s
    @param u Wind vector in U direction (m/s)
    @param v Wind vector in V direction (m/s)
    @return Wind direction (degrees)
    """
    return ((180 / math.pi * math.atan2(-u, -v)) + 360) % 360


def kelvin_to_celcius(t):
    """!
    Convert temperature in Kelvin to Celcius
    @param t Temperature (Celcius)
    @return Temperature (Kelvin)
    """
    return t - 273.15

c = cdsapi.Client(verify=False)

def get_year(year = 2008):
    csv_file = os.path.join(DIR, '{:04d}.csv'.format(year))

    def to_datetime(hours):
        return(datetime.datetime(1900, 1, 1) + datetime.timedelta(hours=hours))

    MONTHS = list(map(lambda x: '{:02d}'.format(x), range(1, 13)))
    DAYS = list(map(lambda x: '{:02d}'.format(x), range(1, 32)))
    area = "/".join(map(str, AREA))
    time = list(map(lambda x: '{:02d}:00'.format(x), range(0, 24)))
    file_by_year = os.path.join(DIR, '{year:04d}.nc'.format(year=year))
    # file_temp_rh = os.path.join(DIR, '{year:04d}_temp_rh.nc'.format(year=year))
    file_temp = os.path.join(DIR, '{year:04d}_temp.nc'.format(year=year))
    file_dew = os.path.join(DIR, '{year:04d}_dew.nc'.format(year=year))
    # file_wind = os.path.join(DIR, '{year:04d}_wind.nc'.format(year=year))
    file_wind_u = os.path.join(DIR, '{year:04d}_wind_u.nc'.format(year=year))
    file_wind_v = os.path.join(DIR, '{year:04d}_wind_v.nc'.format(year=year))
    file_prec = os.path.join(DIR, '{year:04d}_prec.nc'.format(year=year))
    if not os.path.exists(file_temp):
        c.retrieve(
            'reanalysis-era5-land',
            {
                'product_type': 'reanalysis',
                'variable': '2m_temperature',
                'year': '{:04d}'.format(year),
                'month': MONTHS,
                'day': DAYS,
                'time': time,
                'area': AREA,
                'format': 'netcdf',
            },
            file_temp)
    if not os.path.exists(file_dew):
        c.retrieve(
            'reanalysis-era5-land',
            {
                'product_type': 'reanalysis',
                'variable': '2m_dewpoint_temperature',
                'year': '{:04d}'.format(year),
                'month': MONTHS,
                'day': DAYS,
                'time': time,
                'area': AREA,
                'format': 'netcdf',
            },
            file_dew)
    if not os.path.exists(file_wind_u):
        c.retrieve(
            'reanalysis-era5-land',
            {
                'product_type': 'reanalysis',
                'variable': '10m_u_component_of_wind',
                'year': '{:04d}'.format(year),
                'month': MONTHS,
                'day': DAYS,
                'time': time,
                'area': AREA,
                'format': 'netcdf',
            },
            file_wind_u)
    if not os.path.exists(file_wind_v):
        c.retrieve(
            'reanalysis-era5-land',
            {
                'product_type': 'reanalysis',
                'variable': '10m_v_component_of_wind',
                'year': '{:04d}'.format(year),
                'month': MONTHS,
                'day': DAYS,
                'time': time,
                'area': AREA,
                'format': 'netcdf',
            },
            file_wind_v)
    if not os.path.exists(file_prec):
        c.retrieve(
            'reanalysis-era5-land',
            {
                'product_type': 'reanalysis',
                'variable': ['total_precipitation'],
                'year': '{:04d}'.format(year),
                'month': MONTHS,
                'day': DAYS,
                'time': time,
                'area': AREA,
                'format': 'netcdf',
            },
            file_prec)

    nc_temp = netCDF4.Dataset(file_temp)
    nc_dew = netCDF4.Dataset(file_dew)
    nc_wind_u = netCDF4.Dataset(file_wind_u)
    nc_wind_v = netCDF4.Dataset(file_wind_v)
    nc_prec = netCDF4.Dataset(file_prec)
    start_date = to_datetime(int(nc_temp.variables['time'][[0]].data[0]))
    df = pd.DataFrame(data=None, columns=["date", "latitude", "longitude", "temp", "rh", "ws", "wd", "prec"])
    df.to_csv(csv_file, index=False)
    for k in ['time', 'latitude', 'longitude']:
        assert(all(nc_temp.variables[k][:].data == nc_dew.variables[k][:].data))
        assert(all(nc_temp.variables[k][:].data == nc_wind_u.variables[k][:].data))
        assert(all(nc_temp.variables[k][:].data == nc_wind_v.variables[k][:].data))
        assert(all(nc_temp.variables[k][:].data == nc_prec.variables[k][:].data))
    #result = []

    for i_time, time in enumerate(nc_temp.variables['time']):
        date = to_datetime(int(time.data))
        print(date)
        by_latitude = []
        for i_lat, latitude in enumerate(nc_temp.variables['latitude']):
            longitudes = nc_temp.variables['longitude'][:].data
            # t = nc_temp_rh.variables['t']
            # r = nc_temp_rh.variables['r']
            # tp = nc_prec.variables['tp']
            # u = nc_wind.variables['u']
            # v = nc_wind.variables['v']
            # temp = t[i_time, i_lat].data * t.scale_factor + t.add_offset
            # rh = r[i_time, i_lat].data * r.scale_factor + r.add_offset
            # prec = tp[i_time, i_lat].data * tp.scale_factor + tp.add_offset
            # wind_u = u[i_time, i_lat].data * u.scale_factor + u.add_offset
            # wind_v = v[i_time, i_lat].data * v.scale_factor + v.add_offset
            def get_var(nc, name):
                v = nc.variables[name]
                d = v[i_time, i_lat].data
                def fix(x):
                    return(np.NaN if x == v.missing_value else x)
                d = np.vectorize(fix)(d)
                return(d)
            temp = get_var(nc_temp, 't2m')
            dew = get_var(nc_dew, 'd2m')
            prec = get_var(nc_prec, 'tp')
            u = get_var(nc_wind_u, 'u10')
            v = get_var(nc_wind_v, 'v10')
            df = pd.DataFrame({'longitude': longitudes,
                               'temp': np.vectorize(kelvin_to_celcius)(temp),
                               'rh': np.vectorize(calc_rh)(np.vectorize(kelvin_to_celcius)(dew),
                                                           np.vectorize(kelvin_to_celcius)(temp)),
                               'ws': np.vectorize(calc_ws)(u, v),
                               'wd': np.vectorize(calc_wd)(u, v),
                               # convert m to mm
                               'prec': prec * 1000})
            df['latitude'] = latitude
            df = df[np.isfinite(df['temp'])]
            df = df[np.isfinite(df['rh'])]
            df = df[np.isfinite(df['ws'])]
            df = df[np.isfinite(df['wd'])]
            by_latitude.append(df)
        df = pd.concat(by_latitude)
        df['time'] = date
        #result.append(df)
        df = df[['time', 'latitude', 'longitude', 'temp', 'rh', 'ws', 'wd', 'prec']]
        df['temp'] = df['temp'].apply(lambda x: '{:0.1f}'.format(x))
        df['ws'] = df['ws'].apply(lambda x: '{:0.1f}'.format(x))
        df['wd'] = df['wd'].apply(lambda x: '{:0.0f}'.format(x))
        df['prec'] = df['prec'].apply(lambda x: '{:0.1f}'.format(x))
        df.to_csv(csv_file, header=False, index=False, mode='a')

for y in range(1980, 2021):
    get_year(y)
