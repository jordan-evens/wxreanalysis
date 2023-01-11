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

stns <- list()
# for (s in list.files(DIR))
# {
#   stn <- substr(s, 1, 3)
#   df <- readJSON(paste0(DIR, s))
#   stns[[stn]] <- df
# }

stations <- as.data.table(read.csv('stations.csv'))


for (i in 1:nrow(stations))
{
  r <- stations[i]
  if ('Spot Check' == r$Analysis)
  {
    df <- readJSON(paste0(DIR, r$station, '.json'))
    stns[[r$station]] <- df
  }
}

COLS_WX <- c('ID', 'LAT', 'LONG', 'DATE', 'TEMP', 'RH', 'WS', 'WD', 'PREC')
COLS_FWI <- c('FFMC', 'DMC', 'DC', 'ISI', 'BUI', 'FWI', 'DSR')
COLS_ALL <- c(COLS_WX, COLS_FWI)

get_AB <- function(filename='C:/nrcan/data/AB_AAF_2000-2017_wxObs.csv') {
  wx_ab <- data.table::fread(filename)
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

calc_all <- function(dir_in, wx) {
  wx_fwi_all <- NULL
  wx_orig_all <- NULL
  wx_recalc_all <- NULL
  for (i in 1:length(MODELS)) {
    model <- MODELS[i]
    DIR <- paste0(dir_in, '/', model, '/')
    for (stn in names(stns)) {
      print(stn)
      wx_model <- stns[[stn]]
      wx_model[, `:=`(temp = as.double(temp), rh = as.double(rh), ws = as.double(ws), wd=as.integer(wd), prec=as.double(prec), latitude=as.double(latitude), longitude=as.double(longitude))]
      wx_stn <- wx_on[ID == stn]
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
          wx_fwi <- wx_fwi[, ..c(list('MODEL'), COLS_ALL)]
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
  return(list(wx_fwi_all, wx_orig_all, wx_recalc_all))
}


wx_ab <- get_AB()
wx_on <- get_ON()

wx_fwi_on <- calc_all('wx', wx_on)

