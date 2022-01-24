
import ssl
from tqdm import tqdm

## So HTTPS transfers work properly
ssl._create_default_https_context = ssl._create_unverified_context

import requests
from urllib3.exceptions import InsecureRequestWarning
import os

# DIR = './merra2'
DIR = r'G:\wxreanalysis\merra2'
if not os.path.exists(DIR):
    os.makedirs(DIR)

# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)
import urllib3
urllib3.disable_warnings(category=InsecureRequestWarning)


# # LATITUDE_MIN = 41
# # LATITUDE_MAX = 84
# # LONGITUDE_MIN = -141
# # LONGITUDE_MAX = -52
# # AREA = [LATITUDE_MAX, LONGITUDE_MIN, LATITUDE_MIN, LONGITUDE_MAX]
# #
# # # https://disc.gsfc.nasa.gov/datasets/M2I1NXLFO_5.12.4/summary
# # # MERRA-2 inst1_2d_lfo_Nx: 2d,1-Hourly,Instantaneous,Single-Level,Assimilation,
# # # Land Surface Forcings V5.12.4 (M2I1NXLFO)
# # #BASE = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/"
# # collection_shortname = 'M2T1NXAER'
# # collection_longname  = 'tavg1_2d_aer_Nx'
# # collection_number = 'MERRA2_400'
# # MERRA2_version = '5.12.4'
# # BASE = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/{}.{}/{}/{:0>2d}'.format(collection_shortname, MERRA2_version, year, month)
# # #DIR = "MERRA2/M2I1NXLFO.5.12.4/"
# # year = 1982
# # month = 1
# # day = 1
# # ymd = "{:04d}{:02d}{:02d}".format(year, month, day)
# # hours = "[0:23]"
# # lat = "[262:348]"
# # lon = "[62:205]"
# # vars = ["QLML", "SPEEDLML", "TLML"]
# # var_string = ",".join(list(map(lambda x:
# #                                "{var:s}{hours:s}{lat:s}{lon:s}".format(var=x,
# #                                                                        hours=hours,
# #                                                                        lat=lat,
# #                                                                        lon=lon),
# #                                vars)))
# # #url = BASE + DIR + "{year:04d}/{month:02d}/MERRA2_100.inst1_2d_lfo_Nx.{ymd:s}.nc4.nc4?{vars:s}time,lat{lat:s},lon{lon:s}".format(year=year, month=month, ymd=ymd, lat=lat, lon=lon, hours=hours, vars=var_string)
# # to_dir = "./merra-2"
# # username = "jordan_evens_nrcan"
# # password = "eQnhX29-kw!Y/yA"
# # from pydap.client import open_url
# # from pydap.cas.urs import setup_session
# # #dataset_url = 'https://disc.gsfc.nasa.gov/datasets/M2I1NXLFO_5.12.4/'
# # session = setup_session(username, password, check_url=BASE)
# # dataset = open_url(BASE, session=session)
# # ----------------------------------
# # Import Python modules
# # ----------------------------------
#
# import warnings
#
# warnings.filterwarnings("ignore")
#
# import xarray as xr
# # import matplotlib.pyplot as plt
# # import cartopy.crs as ccrs
# from calendar import monthrange
# import time
# import platform
# print('python_version() is ', platform.python_version())
# # ---------------------------------
# # Read data
# # ---------------------------------
# # MERRA-2 collection (hourly)
# collection_shortname = 'M2T1NXAER'
# collection_longname = 'tavg1_2d_aer_Nx'
# collection_number = 'MERRA2_400'
# MERRA2_version = '5.12.4'
# year = 2020
#
# # Open dataset
# # Read selected days in the same month and year
# month = 1  # January
# day_beg = 1
# day_end = 31
#
# # Note that collection_number is MERRA2_401 in a few cases, refer to "Records of MERRA-2 Data Reprocessing and Service Changes"
# if year == 2020 and month == 9:
#     collection_number = 'MERRA2_401'
#
# # OPeNDAP URL
#
# url = 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/{}.{}/{}/{:0>2d}'.format(collection_shortname,
#                                                                                         MERRA2_version, year, month)
# files_month = ['{}/{}.{}.{}{:0>2d}{:0>2d}.nc4'.format(url, collection_number, collection_longname, year, month, days)
#                for days in range(day_beg, day_end + 1, 1)]
# # Get the number of files
# len_files_month = len(files_month)
#
# # Print
# print("{} files to be opened:".format(len_files_month))
# print("files_month", files_month)
#
# print(" ")
# print("Opening...(It may take ~ 5 minutes to open 1-month data)")
# print(" ")
#
# # Read dataset URLs
# #ds = xr.open_mfdataset(files_month)
#
# # View metadata (function like ncdump -c)
# #ds
# url = files_month[0]
# url = "https://acdisc.gesdisc.eosdis.nasa.gov/data//Aqua_AIRS_Level3/AIRX3STD.006/2006/AIRS.2006.12.31.L3.RetStd001.v6.0.9.0.G13155192744.hdf"
# username = "jordan_evens_nrcan"
# password = "eQnhX29-kw!Y/yA"
# from pydap.client import open_url
# from pydap.cas.urs import setup_session
# # session = setup_session(username, password, check_url=url)
# # dataset = open_url(url, session=session)

username = "jordan_evens_nrcan"
password = "eQnhX29-kw!Y/yA"

url_main = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/{year:04d}/{month:02d}/MERRA2_{version:03d}.inst1_2d_asm_Nx.{year:04d}{month:02d}{day:02d}.nc4"
url_flux = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/{year:04d}/{month:02d}/MERRA2_{version:03d}.tavg1_2d_flx_Nx.{year:04d}{month:02d}{day:02d}.nc4"
import os
import requests

# overriding requests.Session.rebuild_auth to mantain headers when redirected

class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = 'urs.earthdata.nasa.gov'
    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)
    # Overrides from the library to keep headers when redirected to or from
    # the NASA auth host.
    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url
        if 'Authorization' in headers:
            original_parsed = requests.utils.urlparse(response.request.url)
            redirect_parsed = requests.utils.urlparse(url)
            if ((original_parsed.hostname != redirect_parsed.hostname)
                    and redirect_parsed.hostname != self.AUTH_HOST
                    and original_parsed.hostname != self.AUTH_HOST):
                del headers['Authorization']
        return

def download(url, filename=None):
    if filename is None:
        filename = os.path.join(DIR, os.path.basename(url))
    if not os.path.exists(filename):
        session = SessionWithHeaderRedirection(username, password)
        try:
            print(filename)
            # submit the request using the session
            response = session.get(url, stream=True, verify=False)
            # print(response.status_code)
            # raise an exception in case of http errors
            response.raise_for_status()
            size = int(response.headers.get('content-length', 0))
            # save the file
            progress_bar = tqdm(total=size, unit='iB', unit_scale=True)
            with open(filename, 'wb') as fd:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    progress_bar.update(len(chunk))
                    fd.write(chunk)
            progress_bar.close()
        except requests.exceptions.HTTPError as e:
            # handle any errors here
            print(e)
    return(filename)

def list_vars(nc):
    return(list(map(lambda x: x + ": " + nc[x].long_name, nc.variables)))

import netCDF4
import datetime

def get_date(date):
    print(date)
    year = date.year
    month = date.month
    day = date.day
    version = 100 if year < 1992 else 200 if year < 2001 else 300 if year < 2011 else 400
    # ['lon: longitude',
    #  'lat: latitude',
    #  'time: time',
    #  'DISPH: zero_plane_displacement_height',
    #  'PS: surface_pressure',
    #  'QV10M: 10-meter_specific_humidity',
    #  'QV2M: 2-meter_specific_humidity',
    #  'SLP: sea_level_pressure',
    #  'T10M: 10-meter_air_temperature',
    #  'T2M: 2-meter_air_temperature',
    #  'TO3: total_column_ozone',
    #  'TOX: total_column_odd_oxygen',
    #  'TQI: total_precipitable_ice_water',
    #  'TQL: total_precipitable_liquid_water',
    #  'TQV: total_precipitable_water_vapor',
    #  'TROPPB: tropopause_pressure_based_on_blended_estimate',
    #  'TROPPT: tropopause_pressure_based_on_thermal_estimate',
    #  'TROPPV: tropopause_pressure_based_on_EPV_estimate',
    #  'TROPQ: tropopause_specific_humidity_using_blended_TROPP_estimate',
    #  'TROPT: tropopause_temperature_using_blended_TROPP_estimate',
    #  'TS: surface_skin_temperature',
    #  'U10M: 10-meter_eastward_wind',
    #  'U2M: 2-meter_eastward_wind',
    #  'U50M: eastward_wind_at_50_meters',
    #  'V10M: 10-meter_northward_wind',
    #  'V2M: 2-meter_northward_wind',
    #  'V50M: northward_wind_at_50_meters']
    try:
        download(url_main.format(year=year, month=month, day=day, version=version))
    except KeyboardInterrupt as e:
        raise e
    except:
        pass
    # ['lon: longitude',
    #  'lat: latitude',
    #  'time: time',
    #  'BSTAR: surface_bouyancy_scale',
    #  'CDH: surface_exchange_coefficient_for_heat',
    #  'CDM: surface_exchange_coefficient_for_momentum',
    #  'CDQ: surface_exchange_coefficient_for_moisture',
    #  'CN: surface_neutral_drag_coefficient',
    #  'DISPH: zero_plane_displacement_height',
    #  'EFLUX: total_latent_energy_flux',
    #  'EVAP: evaporation_from_turbulence',
    #  'FRCAN: areal_fraction_of_anvil_showers',
    #  'FRCCN: areal_fraction_of_convective_showers',
    #  'FRCLS: areal_fraction_of_nonanvil_large_scale_showers',
    #  'FRSEAICE: ice_covered_fraction_of_tile',
    #  'GHTSKIN: Ground_heating_for_skin_temp',
    #  'HFLUX: sensible_heat_flux_from_turbulence',
    #  'HLML: surface_layer_height',
    #  'NIRDF: surface_downwelling_nearinfrared_diffuse_flux',
    #  'NIRDR: surface_downwelling_nearinfrared_beam_flux',
    #  'PBLH: planetary_boundary_layer_height',
    #  'PGENTOT: Total_column_production_of_precipitation',
    #  'PRECANV: anvil_precipitation',
    #  'PRECCON: convective_precipitation',
    #  'PRECLSC: nonanvil_large_scale_precipitation',
    #  'PRECSNO: snowfall',
    #  'PRECTOT: total_precipitation',
    #  'PRECTOTCORR: total_precipitation',
    #  'PREVTOT: Total_column_re-evap/subl_of_precipitation',
    #  'QLML: surface_specific_humidity',
    #  'QSH: effective_surface_specific_humidity',
    #  'QSTAR: surface_moisture_scale',
    #  'RHOA: air_density_at_surface',
    #  'RISFC: surface_bulk_richardson_number',
    #  'SPEED: surface_wind_speed',
    #  'SPEEDMAX: surface_wind_speed',
    #  'TAUGWX: surface_eastward_gravity_wave_stress',
    #  'TAUGWY: surface_northward_gravity_wave_stress',
    #  'TAUX: eastward_surface_stress',
    #  'TAUY: northward_surface_stress',
    #  'TCZPBL: transcom_planetary_boundary_layer_height',
    #  'TLML: surface_air_temperature',
    #  'TSH: effective_surface_skin_temperature',
    #  'TSTAR: surface_temperature_scale',
    #  'ULML: surface_eastward_wind',
    #  'USTAR: surface_velocity_scale',
    #  'VLML: surface_northward_wind',
    #  'Z0H: surface_roughness_for_heat',
    #  'Z0M: surface_roughness']
    try:
        download(url_flux.format(year=year, month=month, day=day, version=version))
    except KeyboardInterrupt as e:
        raise e
    except:
        pass

d = datetime.date(1980, 1, 1)
while (d < datetime.date.today()):
    get_date(d)
    d = d + datetime.timedelta(days=1)
