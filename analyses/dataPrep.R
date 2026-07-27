###################Started by Mao####################
###################July 23,2026######################

library("dplR")
library("dplyr")
setwd("C:/PhD/Project/PhD_thesis/mast_growth")

###Read in the ring width data
AB08THSE <- read.rwl("data/measurement/crossDated/AB08_TSHE_dated.rwl", format = "auto")

###Read in the DBH data
dbh <- read.csv("data/dbhMORA.csv",header = TRUE)

mora_stand <- c("AB08", "AV06", "TO04", "TA01", "AO03", "AG05", "AE10", "AM16")
dbh <- dbh[dbh$STANDID %in% mora_stand, ]

latest_dbh <- do.call(rbind, lapply(split(dbh, list(dbh$STANDID, dbh$TAG)), function(x) {
  x[which.max(x$YEAR), ]
}))

