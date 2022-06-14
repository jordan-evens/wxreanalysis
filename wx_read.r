library(rjson)
library(data.table)
library(lubridate)

DIR <- 'wx/'
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

for (t in unique(stations$Analysis))
{
  print(t)
  if (t != 'Spot Check')
  {
    by_type <- stations[stations$Analysis == t]
    # result <- NULL
    # for (stn in unique(by_type$station))
    # {
    #   result <- rbind(result, stns[[stn]])
    # }
    # result$distance <- NULL
    # result <- unique(result)
    # stns[[t]] <- result
    df <- readJSON(paste0(DIR, t, '.json'))
    stns[[t]] <- df
  }
}

for (i in 1:nrow(stations))
{
  r <- stations[i]
  if ('Spot Check' == r$Analysis)
  {
    df <- readJSON(paste0(DIR, r$station, '.json'))
    stns[[r$station]] <- df
  }
}