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
import json
import pytz

from wx_era5 import *

df = pd.read_csv('stations.csv', index_col='FID')
out_dir = ensure_dir('wx')
for stn in df.itertuples():
    print(stn)
    start = datetime.datetime(stn.minYR, 1, 1, tzinfo=pytz.timezone('GMT'))
    end = datetime.datetime(stn.maxYr + 1, 1, 1, tzinfo=pytz.timezone('GMT'))
    hours = (end - start).days * 24 - 1
    # HACK: daylight savings or leap years are making a mess of things?
    while((start + datetime.timedelta(hours=hours + 1)) <= end):
        hours = hours + 1
    if stn.Analysis == 'Spot Check':
        out_file = os.path.join(out_dir, '{}.json'.format(stn.station))
        if not os.path.exists(out_file):
            wx = getWeather(stn.latitude, stn.longitude, 'ERA5', 'GMT', N=5, start=start, hours=hours)
            with open(out_file, 'w') as f:
                json.dump(json.loads(wx), f, indent=2)

for a in df['Analysis'].unique():
    if a != 'Spot Check':
        print(a)
        cur = df[df['Analysis'] == a]
        start = datetime.datetime(max(cur.minYR), 1, 1, tzinfo=pytz.timezone('GMT'))
        end = datetime.datetime(min(cur.maxYr) + 1, 1, 1, tzinfo=pytz.timezone('GMT'))
        hours = (end - start).days * 24 - 1
        # HACK: daylight savings or leap years are making a mess of things?
        while ((start + datetime.timedelta(hours=hours + 1)) <= end):
            hours = hours + 1
        bounds = {
            'lat_max': max(cur.latitude) + 0.5,
            'lon_min': min(cur.longitude) - 0.5,
            'lon_max': max(cur.longitude) + 0.5,
            'lat_min': min(cur.latitude) - 0.5,
        }
        out_file = os.path.join(out_dir, '{}.json'.format(a))
        if not os.path.exists(out_file):
            wx = getWeatherByBounds(bounds, 'ERA5', 'GMT', start=start, hours=hours)
            with open(out_file, 'w') as f:
                json.dump(json.loads(wx), f, indent=2)
