library(dplyr)
library(ggplot2)
library(reshape2)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/')
filePath <- "All_Subjects_Pepisode_21May2015_CLEAN.csv"
pepisode <- read.csv(filePath)

# Initialize output and row counter
statsOutput <- data.frame(FrequencyBand = NA, Comparison = NA, PValue = NA)
thisRow <- 1

# get rid of NaN
pepisode <- pepisode %>%
  filter(is.na(Pepisode) == FALSE)

# put frequency bands in order
freqBandOrder <- c("Delta","Theta","Alpha","Beta","Gamma")
pepisode$FrequencyBand <- factor(pepisode$FrequencyBand, levels = freqBandOrder)

# put time points in order
timePointOrder <- c("Pre","Tele","Post")
pepisode$TimePoint <- factor(pepisode$TimePoint, levels = timePointOrder)

for (thisFreq in 1:nlevels(pepisode$FrequencyBand)){
  freqBandData <- pepisode %>%
    select(TimePoint, Pepisode:FrequencyBand) %>%
    filter(FrequencyBand == levels(pepisode$FrequencyBand)[thisFreq]) %>%
    group_by(ElectrodeID, TimePoint) %>%
    summarise(Pepisode = mean(Pepisode))
  
  # cast to wide format
  freqBandDataWide <- dcast(freqBandData, ElectrodeID ~ TimePoint, value.var = "Pepisode")
  
  # Wilcoxon signed-rank tests
  preVsTele  <- wilcox.test(freqBandDataWide$Pre, freqBandDataWide$Tele, paired = TRUE)
  teleVsPost <- wilcox.test(freqBandDataWide$Tele, freqBandDataWide$Post, paired = TRUE)
  
  # Save the results
  statsOutput[thisRow,]     <- c(levels(pepisode$FrequencyBand)[thisFreq], "Pre vs Tele", preVsTele$p.value)
  statsOutput[thisRow + 1,] <- c(levels(pepisode$FrequencyBand)[thisFreq], "Tele vs Post", teleVsPost$p.value)
  
  thisRow <- thisRow + 2
}



