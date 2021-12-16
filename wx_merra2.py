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

# TLML = temp
# QLML = specific humidity
url4 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXASM.5.12.4/1980/01/MERRA2_100.inst1_2d_asm_Nx.19800101.nc4"
url5 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXINT.5.12.4/1980/01/MERRA2_100.inst1_2d_int_Nx.19800101.nc4"
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I1NXLFO.5.12.4/1980/01/MERRA2_100.inst1_2d_lfo_Nx.19800101.nc4"
url6 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NXGAS.5.12.4/1980/01/MERRA2_100.inst3_2d_gas_Nx.19800101.nc4"
url2 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2SDNXSLV.5.12.4/1980/01/MERRA2_100.statD_2d_slv_Nx.19800101.nc4"
url3 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXADG.5.12.4/1980/01/MERRA2_100.tavg1_2d_adg_Nx.19800101.nc4"
url7 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXAER.5.12.4/1980/01/MERRA2_100.tavg1_2d_aer_Nx.19800101.nc4"
url8 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXSLV.5.12.4/1980/01/MERRA2_100.tavg1_2d_slv_Nx.19800101.nc4"
url9 = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T1NXFLX.5.12.4/1980/01/MERRA2_100.tavg1_2d_flx_Nx.19800101.nc4"
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

def download(filename, url):
    if not os.path.exists(filename):
        session = SessionWithHeaderRedirection(username, password)
        try:
            # submit the request using the session
            response = session.get(url, stream=True)
            print(response.status_code)
            # raise an exception in case of http errors
            response.raise_for_status()
            # save the file
            with open(filename, 'wb') as fd:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    fd.write(chunk)
        except requests.exceptions.HTTPError as e:
            # handle any errors here
            print(e)
    return(filename)

def list_vars(nc):
    return(list(map(lambda x: x + ": " + nc[x].long_name, nc.variables)))

import netCDF4
# ['lon: longitude',
#  'lat: latitude',
#  'time: time',
#  'HLML: surface_layer_height',
#  'PS: surface_pressure',
#  'QLML: surface_specific_humidity',
#  'SPEEDLML: surface_wind_speed',
#  'TLML: surface_air_temperature']
nc1 = netCDF4.Dataset(download("test.nc", url))
# ['lon: longitude',
#  'lat: latitude',
#  'time: time',
#  'HOURNORAIN: time-during_an_hour_with_no_precipitation',
#  'T2MMAX: 2-meter_air_temperature',
#  'T2MMEAN: 2-meter_air_temperature',
#  'T2MMIN: 2-meter_air_temperature',
#  'TPRECMAX: total_precipitation']
nc2 = netCDF4.Dataset(download("test2.nc", url2))
# ['lon: longitude',
#  'lat: latitude',
#  'time: time',
#  'BCDP001: Black Carbon Dry Deposition Bin 001',
#  'BCDP002: Black Carbon Dry Deposition Bin 002',
#  'BCEM001: Black Carbon Emission Bin 001',
#  'BCEM002: Black Carbon Emission Bin 002',
#  'BCEMAN: Black Carbon Anthropogenic Emissions',
#  'BCEMBB: Black Carbon Biomass Burning Emissions',
#  'BCEMBF: Black Carbon Biofuel Emissions',
#  'BCHYPHIL: Black Carbon Hydrophobic to Hydrophilic',
#  'BCSD001: Black Carbon Sedimentation Bin 001',
#  'BCSD002: Black Carbon Sedimentation Bin 002',
#  'BCSV001: Black Carbon Convective Scavenging Bin 001',
#  'BCSV002: Black Carbon Convective Scavenging Bin 002',
#  'BCWT001: Black Carbon Wet Deposition Bin 001',
#  'BCWT002: Black Carbon Wet Deposition Bin 002',
#  'DUAERIDX: Dust TOMS UV Aerosol Index',
#  'DUDP001: Dust Dry Deposition Bin 001',
#  'DUDP002: Dust Dry Deposition Bin 002',
#  'DUDP003: Dust Dry Deposition Bin 003',
#  'DUDP004: Dust Dry Deposition Bin 004',
#  'DUDP005: Dust Dry Deposition Bin 005',
#  'DUEM001: Dust Emission Bin 001',
#  'DUEM002: Dust Emission Bin 002',
#  'DUEM003: Dust Emission Bin 003',
#  'DUEM004: Dust Emission Bin 004',
#  'DUEM005: Dust Emission Bin 005',
#  'DUEXTTFM: Dust Extinction AOT [550 nm] - PM 1.0 um',
#  'DUSCATFM: Dust Scattering AOT [550 nm] - PM 1.0 um',
#  'DUSD001: Dust Sedimentation Bin 001',
#  'DUSD002: Dust Sedimentation Bin 002',
#  'DUSD003: Dust Sedimentation Bin 003',
#  'DUSD004: Dust Sedimentation Bin 004',
#  'DUSD005: Dust Sedimentation Bin 005',
#  'DUSV001: Dust Convective Scavenging Bin 001',
#  'DUSV002: Dust Convective Scavenging Bin 002',
#  'DUSV003: Dust Convective Scavenging Bin 003',
#  'DUSV004: Dust Convective Scavenging Bin 004',
#  'DUSV005: Dust Convective Scavenging Bin 005',
#  'DUWT001: Dust Wet Deposition Bin 001',
#  'DUWT002: Dust Wet Deposition Bin 002',
#  'DUWT003: Dust Wet Deposition Bin 003',
#  'DUWT004: Dust Wet Deposition Bin 004',
#  'DUWT005: Dust Wet Deposition Bin 005',
#  'OCDP001: Organic Carbon Dry Deposition Bin 001 __ENSEMBLE__',
#  'OCDP002: Organic Carbon Dry Deposition Bin 002 __ENSEMBLE__',
#  'OCEM001: Organic Carbon Emission Bin 001 __ENSEMBLE__',
#  'OCEM002: Organic Carbon Emission Bin 002 __ENSEMBLE__',
#  'OCEMAN: Organic Carbon Anthropogenic Emissions __ENSEMBLE__',
#  'OCEMBB: Organic Carbon Biomass Burning Emissions __ENSEMBLE__',
#  'OCEMBF: Organic Carbon Biofuel Emissions __ENSEMBLE__',
#  'OCEMBG: Organic Carbon Biogenic Emissions __ENSEMBLE__',
#  'OCHYPHIL: Organic Carbon Hydrophobic to Hydrophilic __ENSEMBLE__',
#  'OCSD001: Organic Carbon Sedimentation Bin 001 __ENSEMBLE__',
#  'OCSD002: Organic Carbon Sedimentation Bin 002 __ENSEMBLE__',
#  'OCSV001: Organic Carbon Convective Scavenging Bin 001 __ENSEMBLE__',
#  'OCSV002: Organic Carbon Convective Scavenging Bin 002 __ENSEMBLE__',
#  'OCWT001: Organic Carbon Wet Deposition Bin 001 __ENSEMBLE__',
#  'OCWT002: Organic Carbon Wet Deposition Bin 002 __ENSEMBLE__',
#  'SO2EMAN: SO2 Anthropogenic Emissions __ENSEMBLE__',
#  'SO2EMBB: SO2 Biomass Burning Emissions __ENSEMBLE__',
#  'SO2EMVE: SO2 Volcanic (explosive) Emissions __ENSEMBLE__',
#  'SO2EMVN: SO2 Volcanic (non-explosive) Emissions __ENSEMBLE__',
#  'SO4EMAN: SO4 Anthropogenic Emissions __ENSEMBLE__',
#  'SSAERIDX: Sea Salt TOMS UV Aerosol Index',
#  'SSDP001: Sea Salt Dry Deposition Bin 001',
#  'SSDP002: Sea Salt Dry Deposition Bin 002',
#  'SSDP003: Sea Salt Dry Deposition Bin 003',
#  'SSDP004: Sea Salt Dry Deposition Bin 004',
#  'SSDP005: Sea Salt Dry Deposition Bin 005',
#  'SSEM001: Sea Salt Emission Bin 001',
#  'SSEM002: Sea Salt Emission Bin 002',
#  'SSEM003: Sea Salt Emission Bin 003',
#  'SSEM004: Sea Salt Emission Bin 004',
#  'SSEM005: Sea Salt Emission Bin 005',
#  'SSEXTTFM: Sea Salt Extinction AOT [550 nm] - PM 1.0 um',
#  'SSSCATFM: Sea Salt Scattering AOT [550 nm] - PM 1.0 um',
#  'SSSD001: Sea Salt Sedimentation Bin 001',
#  'SSSD002: Sea Salt Sedimentation Bin 002',
#  'SSSD003: Sea Salt Sedimentation Bin 003',
#  'SSSD004: Sea Salt Sedimentation Bin 004',
#  'SSSD005: Sea Salt Sedimentation Bin 005',
#  'SSSV001: Sea Salt Convective Scavenging Bin 001',
#  'SSSV002: Sea Salt Convective Scavenging Bin 002',
#  'SSSV003: Sea Salt Convective Scavenging Bin 003',
#  'SSSV004: Sea Salt Convective Scavenging Bin 004',
#  'SSSV005: Sea Salt Convective Scavenging Bin 005',
#  'SSWT001: Sea Salt Wet Deposition Bin 001',
#  'SSWT002: Sea Salt Wet Deposition Bin 002',
#  'SSWT003: Sea Salt Wet Deposition Bin 003',
#  'SSWT004: Sea Salt Wet Deposition Bin 004',
#  'SSWT005: Sea Salt Wet Deposition Bin 005',
#  'SUDP001: Sulfate Dry Deposition Bin 001 __ENSEMBLE__',
#  'SUDP002: Sulfate Dry Deposition Bin 002 __ENSEMBLE__',
#  'SUDP003: Sulfate Dry Deposition Bin 003 __ENSEMBLE__',
#  'SUDP004: Sulfate Dry Deposition Bin 004 __ENSEMBLE__',
#  'SUEM001: Sulfate Emission Bin 001 __ENSEMBLE__',
#  'SUEM002: Sulfate Emission Bin 002 __ENSEMBLE__',
#  'SUEM003: Sulfate Emission Bin 003 __ENSEMBLE__',
#  'SUEM004: Sulfate Emission Bin 004 __ENSEMBLE__',
#  'SUPMSA: MSA Prod from DMS Oxidation [column] __ENSEMBLE__',
#  'SUPSO2: SO2 Prod from DMS Oxidation [column] __ENSEMBLE__',
#  'SUPSO4AQ: SO4 Prod from Aqueous SO2 Oxidation [column] __ENSEMBLE__',
#  'SUPSO4G: SO4 Prod from Gaseous SO2 Oxidation [column] __ENSEMBLE__',
#  'SUPSO4WT: SO4 Prod from Aqueous SO2 Oxidation (wet dep) [column] __ENSEMBLE__',
#  'SUSD001: Sulfate Settling Bin 001 __ENSEMBLE__',
#  'SUSD002: Sulfate Settling Bin 002 __ENSEMBLE__',
#  'SUSD003: Sulfate Settling Bin 003 __ENSEMBLE__',
#  'SUSD004: Sulfate Settling Bin 004 __ENSEMBLE__',
#  'SUSV001: Sulfate Convective Scavenging Bin 001 __ENSEMBLE__',
#  'SUSV002: Sulfate Convective Scavenging Bin 002 __ENSEMBLE__',
#  'SUSV003: Sulfate Convective Scavenging Bin 003 __ENSEMBLE__',
#  'SUSV004: Sulfate Convective Scavenging Bin 004 __ENSEMBLE__',
#  'SUWT001: Sulfate Wet Deposition Bin 001 __ENSEMBLE__',
#  'SUWT002: Sulfate Wet Deposition Bin 002 __ENSEMBLE__',
#  'SUWT003: Sulfate Wet Deposition Bin 003 __ENSEMBLE__',
#  'SUWT004: Sulfate Wet Deposition Bin 004 __ENSEMBLE__']
nc3 = netCDF4.Dataset(download("test3.nc", url3))
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
nc4 = netCDF4.Dataset(download("test4.nc", url4))
# ['lon: longitude',
#  'lat: latitude',
#  'time: time',
#  'CPT: vertically_integrated_enthalpy',
#  'KE: vertically_integrated_kinetic_energy',
#  'MASS: atmospheric_mass',
#  'THV: vertically_integrated_virtual_potential_temperature',
#  'TOX: total_column_odd_oxygen',
#  'TQI: total_precipitable_ice_water',
#  'TQL: total_precipitable_liquid_water',
#  'TQV: total_precipitable_water_vapor']
nc5 = netCDF4.Dataset(download("test5.nc", url5))
# ['lon: longitude',
#  'lat: latitude',
#  'time: time',
#  'AODANA: Aerosol Optical Depth Analysis',
#  'AODINC: Aerosol Optical Depth Analysis Increment']
nc6 = netCDF4.Dataset(download("test6.nc", url6))
# ['lon: longitude',
#  'lat: latitude',
#  'time: time',
#  'BCANGSTR: Black Carbon Angstrom parameter [470-870 nm]',
#  'BCCMASS: Black Carbon Column Mass Density',
#  'BCEXTTAU: Black Carbon Extinction AOT [550 nm]',
#  'BCFLUXU: Black Carbon column u-wind mass flux',
#  'BCFLUXV: Black Carbon column v-wind mass flux',
#  'BCSCATAU: Black Carbon Scattering AOT [550 nm]',
#  'BCSMASS: Black Carbon Surface Mass Concentration',
#  'DMSCMASS: DMS Column Mass Density __ENSEMBLE__',
#  'DMSSMASS: DMS Surface Mass Concentration __ENSEMBLE__',
#  'DUANGSTR: Dust Angstrom parameter [470-870 nm]',
#  'DUCMASS: Dust Column Mass Density',
#  'DUCMASS25: Dust Column Mass Density - PM 2.5',
#  'DUEXTT25: Dust Extinction AOT [550 nm] - PM 2.5',
#  'DUEXTTAU: Dust Extinction AOT [550 nm]',
#  'DUFLUXU: Dust column u-wind mass flux',
#  'DUFLUXV: Dust column v-wind mass flux',
#  'DUSCAT25: Dust Scattering AOT [550 nm] - PM 2.5',
#  'DUSCATAU: Dust Scattering AOT [550 nm]',
#  'DUSMASS: Dust Surface Mass Concentration',
#  'DUSMASS25: Dust Surface Mass Concentration - PM 2.5',
#  'OCANGSTR: Organic Carbon Angstrom parameter [470-870 nm] __ENSEMBLE__',
#  'OCCMASS: Organic Carbon Column Mass Density __ENSEMBLE__',
#  'OCEXTTAU: Organic Carbon Extinction AOT [550 nm] __ENSEMBLE__',
#  'OCFLUXU: Organic Carbon column u-wind mass flux __ENSEMBLE__',
#  'OCFLUXV: Organic Carbon column v-wind mass flux __ENSEMBLE__',
#  'OCSCATAU: Organic Carbon Scattering AOT [550 nm] __ENSEMBLE__',
#  'OCSMASS: Organic Carbon Surface Mass Concentration __ENSEMBLE__',
#  'SO2CMASS: SO2 Column Mass Density __ENSEMBLE__',
#  'SO2SMASS: SO2 Surface Mass Concentration __ENSEMBLE__',
#  'SO4CMASS: SO4 Column Mass Density __ENSEMBLE__',
#  'SO4SMASS: SO4 Surface Mass Concentration __ENSEMBLE__',
#  'SSANGSTR: Sea Salt Angstrom parameter [470-870 nm]',
#  'SSCMASS: Sea Salt Column Mass Density',
#  'SSCMASS25: Sea Salt Column Mass Density - PM 2.5',
#  'SSEXTT25: Sea Salt Extinction AOT [550 nm] - PM 2.5',
#  'SSEXTTAU: Sea Salt Extinction AOT [550 nm]',
#  'SSFLUXU: Sea Salt column u-wind mass flux',
#  'SSFLUXV: Sea Salt column v-wind mass flux',
#  'SSSCAT25: Sea Salt Scattering AOT [550 nm] - PM 2.5',
#  'SSSCATAU: Sea Salt Scattering AOT [550 nm]',
#  'SSSMASS: Sea Salt Surface Mass Concentration',
#  'SSSMASS25: Sea Salt Surface Mass Concentration - PM 2.5',
#  'SUANGSTR: SO4 Angstrom parameter [470-870 nm] __ENSEMBLE__',
#  'SUEXTTAU: SO4 Extinction AOT [550 nm] __ENSEMBLE__',
#  'SUFLUXU: SO4 column u-wind mass flux __ENSEMBLE__',
#  'SUFLUXV: SO4 column v-wind mass flux __ENSEMBLE__',
#  'SUSCATAU: SO4 Scattering AOT [550 nm] __ENSEMBLE__',
#  'TOTANGSTR: Total Aerosol Angstrom parameter [470-870 nm]',
#  'TOTEXTTAU: Total Aerosol Extinction AOT [550 nm]',
#  'TOTSCATAU: Total Aerosol Scattering AOT [550 nm]']
nc7 = netCDF4.Dataset(download("test7.nc", url7))
# ['lon: longitude',
#  'lat: latitude',
#  'time: time',
#  'CLDPRS: cloud_top_pressure',
#  'CLDTMP: cloud_top_temperature',
#  'DISPH: zero_plane_displacement_height',
#  'H1000: height_at_1000_mb',
#  'H250: height_at_250_hPa',
#  'H500: height_at_500_hPa',
#  'H850: height_at_850_hPa',
#  'OMEGA500: omega_at_500_hPa',
#  'PBLTOP: pbltop_pressure',
#  'PS: surface_pressure',
#  'Q250: specific_humidity_at_250_hPa',
#  'Q500: specific_humidity_at_500_hPa',
#  'Q850: specific_humidity_at_850_hPa',
#  'QV10M: 10-meter_specific_humidity',
#  'QV2M: 2-meter_specific_humidity',
#  'SLP: sea_level_pressure',
#  'T10M: 10-meter_air_temperature',
#  'T250: air_temperature_at_250_hPa',
#  'T2M: 2-meter_air_temperature',
#  'T2MDEW: dew_point_temperature_at_2_m',
#  'T2MWET: wet_bulb_temperature_at_2_m',
#  'T500: air_temperature_at_500_hPa',
#  'T850: air_temperature_at_850_hPa',
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
#  'U250: eastward_wind_at_250_hPa',
#  'U2M: 2-meter_eastward_wind',
#  'U500: eastward_wind_at_500_hPa',
#  'U50M: eastward_wind_at_50_meters',
#  'U850: eastward_wind_at_850_hPa',
#  'V10M: 10-meter_northward_wind',
#  'V250: northward_wind_at_250_hPa',
#  'V2M: 2-meter_northward_wind',
#  'V500: northward_wind_at_500_hPa',
#  'V50M: northward_wind_at_50_meters',
#  'V850: northward_wind_at_850_hPa',
#  'ZLCL: lifting_condensation_level']
nc8 = netCDF4.Dataset(download("test8.nc", url8))
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
nc9 = netCDF4.Dataset(download("test9.nc", url9))
