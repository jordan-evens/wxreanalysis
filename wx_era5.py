#!/usr/bin/env python

import cdsapi
import netCDF4
import datetime
import os
import pandas as pd
import numpy as np
import math
import json
import pytz
from util import calc_rh
from util import calc_wd
from util import calc_ws
from util import ensure_dir
from util import kelvin_to_celcius

LATITUDE_MIN = 41
LATITUDE_MAX = 84
LONGITUDE_MIN = -141
LONGITUDE_MAX = -52
AREA = [LATITUDE_MAX, LONGITUDE_MIN, LATITUDE_MIN, LONGITUDE_MAX]
# NOTE: need to get a client key as per https://cds.climate.copernicus.eu/api-how-to
DATE_MIN = datetime.datetime(1980, 1, 1, tzinfo=pytz.timezone("GMT"))
DATE_MAX = datetime.datetime(2021, 12, 31, 23, tzinfo=pytz.timezone("GMT"))
MODEL = 'ERA5'

DIR = ensure_dir('G:/wxreanalysis/era5')
# DIR = ensure_dir('era5land')

c = cdsapi.Client(verify=False)

def get_year(year = 2008, convert=True):
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
            'reanalysis-era5-single-levels',
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
            'reanalysis-era5-single-levels',
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
            'reanalysis-era5-single-levels',
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
            'reanalysis-era5-single-levels',
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
            'reanalysis-era5-single-levels',
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
    if not convert:
        return
    nc_temp = netCDF4.Dataset(file_temp)
    nc_dew = netCDF4.Dataset(file_dew)
    nc_wind_u = netCDF4.Dataset(file_wind_u)
    nc_wind_v = netCDF4.Dataset(file_wind_v)
    nc_prec = netCDF4.Dataset(file_prec)
    start_date = to_datetime(int(nc_temp.variables['time'][[0]].data[0]))
    df = pd.DataFrame(data=None, columns=["date", "latitude", "longitude", "temp", "rh", "ws", "wd", "prec"])
    csv_file = os.path.join(DIR, '{:04d}.csv'.format(year))
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
                # def fix(x):
                #     return(np.NaN if x == v.missing_value else x)
                # d = np.vectorize(fix)(d)
                d[d == v.missing_value] = np.NaN
                return(d)
            temp = get_var(nc_temp, 't2m')
            dew = get_var(nc_dew, 'd2m')
            prec = get_var(nc_prec, 'tp')
            u = get_var(nc_wind_u, 'u10')
            v = get_var(nc_wind_v, 'v10')
            df = pd.DataFrame({'longitude': longitudes,
                               'temp': kelvin_to_celcius(temp),
                               'rh': calc_rh(kelvin_to_celcius(dew),
                                                           kelvin_to_celcius(temp)),
                               'ws': calc_ws(u, v),
                               'wd': calc_wd(u, v),
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
        df['p'] = df['prec'] - df['prec'].shift(-1)
        #result.append(df)
        df = df[['time', 'latitude', 'longitude', 'temp', 'rh', 'ws', 'wd', 'prec']]
        df['temp'] = df['temp'].apply(lambda x: '{:0.1f}'.format(x))
        df['rh'] = df['rh'].apply(lambda x: '{:0.0f}'.format(x))
        df['ws'] = df['ws'].apply(lambda x: '{:0.1f}'.format(x))
        df['wd'] = df['wd'].apply(lambda x: '{:0.0f}'.format(x))
        df['prec'] = df['prec'].apply(lambda x: '{:0.1f}'.format(x))
        df.to_csv(csv_file, header=False, index=False, mode='a')

# for y in range(DATE_MIN.year, DATE_MAX.year):
#     get_year(y, convert=False)
#     # get_year(y)

# get available points
def getPoints(year):
    file_temp = os.path.join(DIR, '{year:04d}_temp.nc'.format(year=year))
    nc_temp = netCDF4.Dataset(file_temp)
    # bounds would just be halfway to the next point in each direction?
    latitudes = pd.DataFrame(data=nc_temp.variables['latitude'][:].data, columns=['latitude'])
    latitudes['i_lat'] = latitudes.index
    longitudes = pd.DataFrame(data=nc_temp.variables['longitude'][:].data, columns=['longitude'])
    longitudes['i_lon'] = longitudes.index
    latitudes['y_min'] = (latitudes['latitude'] + latitudes['latitude'].shift(-1)) / 2
    latitudes['y_min'] = latitudes['y_min'].apply(lambda x: x if not np.isnan(x) else latitudes['latitude'].min())
    latitudes['y_max'] = (latitudes['latitude'] + latitudes['latitude'].shift(1)) / 2
    latitudes['y_max'] = latitudes['y_max'].apply(lambda x: x if not np.isnan(x) else latitudes['latitude'].max())
    longitudes['x_min'] = (longitudes['longitude'] + longitudes['longitude'].shift(1)) / 2
    longitudes['x_min'] = longitudes['x_min'].apply(lambda x: x if not np.isnan(x) else longitudes['longitude'].min())
    longitudes['x_max'] = (longitudes['longitude'] + longitudes['longitude'].shift(-1)) / 2
    longitudes['x_max'] = longitudes['x_max'].apply(lambda x: x if not np.isnan(x) else longitudes['longitude'].max())
    points = latitudes.merge(longitudes, how='cross')
    return points

def getModels(latitude, longitude, start=datetime.datetime.now(), format='json'):
    result = {'models': [MODEL]}
    if latitude < LATITUDE_MIN or latitude > LATITUDE_MAX:
        result['models'] = []
    if longitude < LONGITUDE_MIN or longitude > LONGITUDE_MAX:
        result['models'] = []
    if start < DATE_MIN or start > DATE_MAX:
        result['models'] = []
    if format != 'json':
        raise NotImplementedError("Only json is supported")
    return json.dumps(result)


def findByDistance(latitude, longitude, year=DATE_MIN.year):
    points = getPoints(year)
    # FIX: this is wrong because the earth isn't flat
    points['distance'] = points.apply(lambda x: math.sqrt(math.pow(x['latitude'] - latitude, 2) + math.pow(x['longitude'] - longitude, 2)), axis=1)
    points = points.sort_values(['distance'])
    return points


def getBounds(latitude, longitude, model, format='json'):
    if latitude < LATITUDE_MIN or latitude > LATITUDE_MAX:
        pass
    if longitude < LONGITUDE_MIN or longitude > LONGITUDE_MAX:
        pass
    # if start < DATE_MIN or start > DATE_MAX:
    #     pass
    if model != MODEL:
        raise ValueError("Unsupported model")
    if format != 'json':
        raise NotImplementedError("Only json is supported")
    points = findByDistance(latitude, longitude)
    closest = points.iloc[0]
    return json.dumps({
        'center': [closest.latitude, closest.longitude],
        'bounds': [closest.x_min, closest.y_min, closest.x_max, closest.x_min]
    })

def getStations(latitude, longitude, model, N=10, offset=0, format='json'):
    if latitude < LATITUDE_MIN or latitude > LATITUDE_MAX:
        raise ValueError("Location out of bounds")
    if longitude < LONGITUDE_MIN or longitude > LONGITUDE_MAX:
        raise ValueError("Location out of bounds")
    if model != MODEL:
        raise ValueError("Unsupported model")
    if format != 'json':
        raise NotImplementedError("Only json is supported")
    points = findByDistance(latitude, longitude)
    points = points.iloc[offset:N + offset]
    result = {}
    for p in points.itertuples():
        result[p.Index] = {
            'location': [p.latitude, p.longitude],
            'distance': p.distance
        }
    return json.dumps(result)

def to_datetime(hours):
    return (datetime.datetime(1900, 1, 1, tzinfo=pytz.timezone("GMT")) + datetime.timedelta(hours=int(hours))).astimezone(pytz.timezone("GMT"))

to_datetime = np.vectorize(to_datetime)

def getWeatherByPoint(i_lat, i_lon, model, zone, ensemble=None, start=datetime.datetime.now(), hours=24, format='json'):
    start = start.astimezone(pytz.timezone(zone))
    end = (start + datetime.timedelta(hours=hours))
    if start < DATE_MIN:
        raise ValueError("Start date out of bounds")
    if end > DATE_MAX:
        raise ValueError("End date out of bounds")
    if model != MODEL:
        raise ValueError("Unsupported model")
    if format != 'json':
        raise NotImplementedError("Only json is supported")
    rows = []
    # do it from day before start so that we can subtract precip
    for year in sorted(set(list(range((start - datetime.timedelta(days=1)).year, end.year)) + [end.year])):
        file_temp = os.path.join(DIR, '{year:04d}_temp.nc'.format(year=year))
        file_dew = os.path.join(DIR, '{year:04d}_dew.nc'.format(year=year))
        file_wind_u = os.path.join(DIR, '{year:04d}_wind_u.nc'.format(year=year))
        file_wind_v = os.path.join(DIR, '{year:04d}_wind_v.nc'.format(year=year))
        file_prec = os.path.join(DIR, '{year:04d}_prec.nc'.format(year=year))
        nc_temp = netCDF4.Dataset(file_temp)
        nc_dew = netCDF4.Dataset(file_dew)
        nc_wind_u = netCDF4.Dataset(file_wind_u)
        nc_wind_v = netCDF4.Dataset(file_wind_v)
        nc_prec = netCDF4.Dataset(file_prec)
        # start_date = to_datetime(int(nc_temp.variables['time'][[0]].data[0]))
        # df = pd.DataFrame(data=None, columns=["date", "latitude", "longitude", "temp", "rh", "ws", "wd", "prec"])
        # csv_file = os.path.join(DIR, '{:04d}.csv'.format(year))
        # df.to_csv(csv_file, index=False)
        for k in ['time', 'latitude', 'longitude']:
            assert (all(nc_temp.variables[k][:].data == nc_dew.variables[k][:].data))
            assert (all(nc_temp.variables[k][:].data == nc_wind_u.variables[k][:].data))
            assert (all(nc_temp.variables[k][:].data == nc_wind_v.variables[k][:].data))
            assert (all(nc_temp.variables[k][:].data == nc_prec.variables[k][:].data))
        time = to_datetime(nc_temp['time'][:])
        temp = nc_temp['t2m'][:, i_lat, i_lon]
        dew = nc_dew['d2m'][:, i_lat, i_lon]
        prec = nc_prec['tp'][:, i_lat, i_lon]
        u = nc_wind_u['u10'][:, i_lat, i_lon]
        v = nc_wind_v['v10'][:, i_lat, i_lon]
        df = pd.DataFrame({'time': time,
                           'temp': kelvin_to_celcius(temp),
                           'rh': calc_rh(kelvin_to_celcius(dew),
                                          kelvin_to_celcius(temp)),
                           'ws': calc_ws(u, v),
                           'wd': calc_wd(u, v),
                           # convert m to mm
                           'prec': prec * 1000})
        rows.append(df)
    df = pd.concat(rows)
    # change precip to not be cumulative
    df['p'] = df['prec'] - df['prec'].shift(1)
    df['prec'] = df.apply(lambda x: x['p'] if x['time'].hour != 1 else x['prec'], axis=1)
    df = df[df.time >= start]
    df = df[df.time < end]
    df = df[['time', 'temp', 'rh', 'ws', 'wd', 'prec']]
    # df['prec_sum'] = df.groupby(pd.to_datetime(df.time.to_numpy()).to_period('d')).prec.cumsum()
    # HACK: can't make them into int directly for some reason
    df['temp'] = df['temp'].apply(lambda x: '{:0.1f}'.format(x)).astype(float)
    df['rh'] = df['rh'].apply(lambda x: '{:0.0f}'.format(x)).astype(float)
    df['ws'] = df['ws'].apply(lambda x: '{:0.1f}'.format(x)).astype(float)
    df['wd'] = df['wd'].apply(lambda x: '{:0.0f}'.format(x)).astype(float)
    df['prec'] = df['prec'].apply(lambda x: '{:0.1f}'.format(x)).astype(float)
    # df['date'] = df['time'].to_period('d')
    result = {
        'location': [float(nc_temp.variables['latitude'][i_lat]), float(nc_temp.variables['longitude'][i_lon])],
        'start': start.astimezone(pytz.timezone(zone)).strftime('%Y-%m-%d %H:%M:00 %Z'),
        'end': end.astimezone(pytz.timezone(zone)).strftime('%Y-%m-%d %H:%M:00 %Z'),
        'hours': hours,
        'temp': list(df['temp'].values),
        'rh': [int(x) for x in df['rh'].values],
        'ws': list(df['ws'].values),
        'wd': [int(x) for x in df['wd'].values],
        'prec': list(df['prec'].values),
    }
    return result

def getWeather(latitude, longitude, model, zone, N=1, ensemble=None, start=datetime.datetime.now(), hours=24, format='json'):
    if latitude < LATITUDE_MIN or latitude > LATITUDE_MAX:
        raise ValueError("Location out of bounds")
    if longitude < LONGITUDE_MIN or longitude > LONGITUDE_MAX:
        raise ValueError("Location out of bounds")
    if start < DATE_MIN:
        raise ValueError("Start date out of bounds")
    start = start.astimezone(pytz.timezone(zone))
    end = (start + datetime.timedelta(hours=hours))
    if end > DATE_MAX:
        raise ValueError("End date out of bounds")
    if model != MODEL:
        raise ValueError("Unsupported model")
    if format != 'json':
        raise NotImplementedError("Only json is supported")
    # select closest  point to requested lat/lon
    points = findByDistance(latitude, longitude)
    results = []
    for i in range(0, N):
        p = points.iloc[i]
        i_lat = int(p.i_lat)
        i_lon = int(p.i_lon)
        wx_stn = getWeatherByPoint(i_lat, i_lon, model, zone, ensemble, start, hours, format)
        wx_stn['distance'] = p.distance
        results.append(wx_stn)
    return json.dumps(results)

def getWeatherByBounds(bounds, model, zone, ensemble=None, start=datetime.datetime.now(), hours=24, format='json'):
    if start < DATE_MIN:
        raise ValueError("Start date out of bounds")
    end = start + datetime.timedelta(hours=hours)
    if end > DATE_MAX:
        raise ValueError("End date out of bounds")
    if model != MODEL:
        raise ValueError("Unsupported model")
    if format != 'json':
        raise NotImplementedError("Only json is supported")
    if bounds is None:
        raise ValueError("No bounds specified")
    if bounds['lat_min'] is None:
        raise ValueError("lat_min not specified")
    if bounds['lat_max'] is None:
        raise ValueError("lat_max not specified")
    if bounds['lon_min'] is None:
        raise ValueError("lon_min not specified")
    if bounds['lon_max'] is None:
        raise ValueError("lon_max not specified")
    # select closest  point to requested lat/lon
    points = findByDistance((bounds['lat_min'] + bounds['lat_max']) / 2.0,
                            (bounds['lon_min'] + bounds['lon_max']) / 2.0)
    results = []
    for i in range(len(points)):
        p = points.iloc[i]
        if (bounds['lat_min'] <= p.latitude <= bounds['lat_max']
                and bounds['lon_min'] <= p.longitude <= bounds['lon_max']):
            i_lat = int(p.i_lat)
            i_lon = int(p.i_lon)
            wx_stn = getWeatherByPoint(i_lat, i_lon, model, zone, ensemble, start, hours, format)
            wx_stn['distance'] = p.distance
            results.append(wx_stn)
    return json.dumps(results)
