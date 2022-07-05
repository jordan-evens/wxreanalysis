#!/usr/bin/env python

import util
import datetime
import os
import pandas as pd
import json
import pytz

# import wx_era5 as wx_model
import wx_era5land as wx_model
# import wx_merra2 as wx_model

df = pd.read_csv('stations.csv', index_col='FID')
out_dir = util.ensure_dir('wx/{}'.format(wx_model.MODEL))
# df = pd.read_csv('ab_ws_selection_final.csv')
# out_dir = util.ensure_dir('wx_AB/{}'.format(wx_model.MODEL))
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

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
            wx = wx_model.getWeather(stn.latitude, stn.longitude, wx_model.MODEL, 'GMT', N=5, start=start, hours=hours)
            with open(out_file, 'w') as f:
                json.dump(json.loads(wx), f, indent=None)

for a in df['Analysis'].unique():
    if a != 'Spot Check':
        print(a)
        cur = df[df['Analysis'] == a]
        start = datetime.datetime(min(cur.minYR), 1, 1, tzinfo=pytz.timezone('GMT'))
        end = datetime.datetime(max(cur.maxYr) + 1, 1, 1, tzinfo=pytz.timezone('GMT'))
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
            wx = wx_model.getWeatherByBounds(bounds, wx_model.MODEL, 'GMT', start=start, hours=hours)
            with open(out_file, 'w') as f:
                json.dump(json.loads(wx), f, indent=None)
