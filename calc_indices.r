library(rjson)
library(data.table)
library(lubridate)
library(cffdrs)

MODELS = c('ERA5', 'ERA5-LAND', 'MERRA2')

readJSON <- function(filename)
{
  j <- fromJSON(file=filename, simplify=FALSE)
  for (i in sequence(length(j)))
  {
    df <- j[[i]]
    tz <- substr(df$start,
                 nchar(df$start) - 2,
                 nchar(df$start))
    df$time <- as.POSIXlt(df$start,tz=tz) +
      (sequence(length(df$prec)) - 1) * 60 * 60
    stopifnot(as.POSIXlt(df$end, tz=tz) == as.POSIXlt(hours(1) + max(df$time), tz=tz))
    df$latitude <- df$location[[1]]
    df$longitude <- df$location[[2]]
    df$location <- NULL
    df$hours <- NULL
    df$start <- NULL
    df$end <- NULL
    j[[i]] <- as.data.table(df)
  }
  result <- rbindlist(j)
  return(result)
}

read_province <- function(dir_in, stations_csv) {
  result <- list()
  stations <- data.table::fread(stations_csv)
  
  for (model in MODELS) {
    print(model)
    stns <- list()
    for (i in 1:nrow(stations)) {
      r <- stations[i]
      if ('Spot Check' == r$Analysis) {
        print(r$station)
        df <- readJSON(paste0(dir_in, '/', model, '/', r$station, '.json'))
        stns[[r$station]] <- df
      }
    }
    result[[model]] <- stns
  }
  return(result)
}

COLS_WX <- c('ID', 'LAT', 'LONG', 'DATE', 'TEMP', 'RH', 'WS', 'WD', 'PREC')
COLS_FWI <- c('FFMC', 'DMC', 'DC', 'ISI', 'BUI', 'FWI', 'DSR')
COLS_ALL <- c(COLS_WX, COLS_FWI)

get_AB <- function(filename='C:/nrcan/data/abwx-new.csv') {
  wx_ab <- data.table::fread(filename)
  wx_ab[, DATE := date(WEATHER_DATE)]
  wx_ab <- wx_ab[, c('STATION_ID', 'LATITUDE', 'LONGITUDE', 'DATE', 'DRY_BULB_TEMPERATURE', 'RELATIVE_HUMIDITY', 'WIND_SPEED_KMH', 'WIND_DIRECTION', 'PRECIPITATION', 'FINE_FUEL_MOISTURE_CODE', 'DUFF_MOISTURE_CODE', 'DROUGHT_CODE', 'BUILD_UP_INDEX', 'INITIAL_SPREAD_INDEX', 'FIRE_WEATHER_INDEX', 'DAILY_SEVERITY_RATING')]
  names(wx_ab) <- COLS_ALL
  wx_ab <- wx_ab[!(is.na(FFMC) | is.na(DMC) | is.na(DC))]
  return(wx_ab)
}


get_ON <- function(filename='C:/nrcan/data/final_weather_archive_1963_2019.csv') {
  wx_on <- data.table::fread(filename)
  setorder(wx_on, 'WSTNID', 'WX_DATE')
  wx_on <- wx_on[, c('WSTNID', 'LATITUDE', 'LONGITUDE', 'WX_DATE', 'TEMP', 'REL_HUM', 'WIND_SPEED', 'WIND_DIR', 'RAIN', 'FFMC', 'DMC', 'DC', 'ISI', 'BUI', 'FWI', 'DSR')]
  names(wx_on) <- COLS_ALL
  return(wx_on)
}

calc_fwi <- function(wx) {
  wx <- copy(wx)
  wx[, `:=`(YR=year(DATE), MON=month(DATE), DAY=day(DATE))]
  wx_fwi <- cffdrs::fwi(input=wx)
  return(wx_fwi[, ..COLS_ALL])
}

calc_all <- function(wx_prov, wx_models) {
  wx_fwi_all <- NULL
  wx_orig_all <- NULL
  wx_recalc_all <- NULL
  for (i in 1:length(MODELS)) {
    model <- MODELS[i]
    stns <- wx_models[[model]]
    for (stn in names(stns)) {
      print(stn)
      stopifnot(stn %in% wx_prov$ID)
      wx_model <- stns[[stn]]
      wx_model[, `:=`(temp = as.double(temp), rh = as.double(rh), ws = as.double(ws), wd=as.integer(wd), prec=as.double(prec), latitude=as.double(latitude), longitude=as.double(longitude))]
      wx_stn <- wx_prov[ID == stn]
      for (yr in unique(year(wx_model$time))) {
        print(yr)
        wx_model_yr <- wx_model[year(time) == yr]
        wx_stn_yr <- wx_stn[year(DATE) == yr]
        wx_model_yr[, for_date := as_date(ifelse(hour(time) > 11, as_date(time) + days(1), as_date(time)))]
        wx_model_yr[, prec := as.double(prec)]
        wx_model_yr[, prec24 := sum(prec), by=list(for_date)]
        wx_model_yr_1200 <- wx_model_yr[hour(time) == 12]
        wx_model_yr_1200[, prec := prec24]
        wx_model_yr_1200 <- wx_model_yr_1200[, -c('for_date', 'prec24')]
        wx_model_yr_1200[, DATE := date(time)]
        # NOTE: startup is values after calculating using actual startup indices, so you'd need to work backwards to find real startup indices
        startup <- wx_stn_yr[min(DATE) == DATE][1]
        shutdown <- wx_stn_yr[max(DATE) == DATE][1]
        wx <- wx_model_yr_1200[DATE >=  startup$DATE & DATE <= shutdown$DATE]
        wx$stn <- stn
        wx <- wx[, c('stn', 'latitude', 'longitude', 'DATE', 'temp', 'rh', 'ws', 'wd', 'prec')]
        names(wx) <- COLS_WX
        # wx_fwi <- cffdrs::fwi(input=wx, init=c(ffmc=startup$FFMC, dmc=startup$DMC, dc=startup$DC))
        # HACK: assume default startup indices for now since startup is wrong anyway
        pts <- unique(wx[, c('LAT', 'LONG')])
        for (j in 1:nrow(pts)) {
          pt <- pts[j,]
          wx_fwi <- calc_fwi(wx[LAT == pt$LAT & LONG == pt$LONG])
          wx_fwi$ID <- paste0(stn, '_', j)
          wx_fwi$MODEL <- model
          cols <- wx_fwi <- cffdrs::fwi(input=wx, init=c(ffmc=startup$FFMC, dmc=startup$DMC, dc=startup$DC))
          wx_fwi <- wx_fwi[, ..cols]
          wx_fwi_all <- rbind(wx_fwi_all, wx_fwi)
        }
        if (1 == i) {
          wx_orig_all <- rbind(wx_orig_all, wx_stn_yr)
          wx_recalc <- calc_fwi(wx_stn_yr[, ..COLS_WX])
          wx_recalc_all <- rbind(wx_recalc_all, wx_recalc)
        }
      }
    }
  }
  return(list(models=wx_fwi_all, orig=wx_orig_all, recalc=wx_recalc_all))
}


wx_ab <- get_AB()
models_ab <- read_province('wx_AB', 'ab_ws_selection_final.csv')
wx_fwi_ab <- calc_all(wx_ab, models_ab)

wx_on <- get_ON()
models_on <- read_province('wx', 'stations.csv')
wx_fwi_on <- calc_all(wx_on, models_on)

