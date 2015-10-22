# Script Name:  1-Load-Data.R
# Author:       Lindsay Vass
# Date:         15 October 2015
# Purpose:      This script will load the csv files output by the Matlab
#               function "m/Run_Calculate_Pepisode.m"

pepisodeFiles <- dir('csv', pattern = '*pepisode*')
powerFiles    <- dir('csv', pattern = '*power*')

pepisodeData <- vector(mode = "list", length = length(pepisodeFiles))
powerData    <- vector(mode = "list", length = length(powerFiles))

for (i in 1:length(pepisodeData)) {
  pepisodeData[[i]] <- read.csv(paste0('csv/', pepisodeFiles[i]))
}
pepisodeData <- data.table::rbindlist(pepisodeData)

for (i in 1:length(powerData)) {
  powerData[[i]] <- read.csv(paste0('csv/', powerFiles[i]))
}
powerData <- data.table::rbindlist(powerData)

save(file = 'Rda/allRawData.Rda', list = c('pepisodeData', 'powerData'))