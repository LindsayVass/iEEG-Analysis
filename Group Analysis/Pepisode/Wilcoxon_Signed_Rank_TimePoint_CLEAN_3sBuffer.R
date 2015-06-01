library(dplyr)
library(ggplot2)
library(reshape2)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/')
load("All_Subjects_Pepisode_CLEAN_3sBuffer2015-05-22.Rda")

# Initialize output and row counter
statsOutput <- data.frame(FrequencyBand = NA, Comparison = NA, PValue = NA)
thisRow <- 1

# get rid of NaN
pepisode <- allPepisodeData %>%
  filter(is.na(Pepisode) == FALSE)

for (thisFreq in 1:nlevels(pepisode$FrequencyBand)){
  freqBandData <- pepisode %>%
    select(TimePoint, Pepisode:FrequencyBand) %>%
    filter(FrequencyBand == levels(pepisode$FrequencyBand)[thisFreq]) %>%
    group_by(ElectrodeID, TimePoint) %>%
    summarise(Pepisode = mean(Pepisode))
  
  # cast to wide format
  freqBandDataWide <- dcast(freqBandData, ElectrodeID ~ TimePoint, value.var = "Pepisode")
  
  # Wilcoxon signed-rank tests
  pre3VsPre2   <- wilcox.test(freqBandDataWide$Pre3, freqBandDataWide$Pre2, paired = TRUE)
  pre2VsPre1   <- wilcox.test(freqBandDataWide$Pre1, freqBandDataWide$Pre2, paired = TRUE)
  pre1VsTele   <- wilcox.test(freqBandDataWide$Pre1, freqBandDataWide$Tele, paired = TRUE)
  teleVsPost1  <- wilcox.test(freqBandDataWide$Tele, freqBandDataWide$Post1, paired = TRUE)
  post1VsPost2 <- wilcox.test(freqBandDataWide$Post1, freqBandDataWide$Post2, paired = TRUE)
  post2VsPost3 <- wilcox.test(freqBandDataWide$Post2, freqBandDataWide$Post3, paired = TRUE)
  
  # Save the results
  statsOutput[thisRow,]     <- c(levels(pepisode$FrequencyBand)[thisFreq], "Pre3 vs Pre2", pre3VsPre2$p.value)
  statsOutput[thisRow + 1,] <- c(levels(pepisode$FrequencyBand)[thisFreq], "Pre2 vs Pre1", pre2VsPre1$p.value)
  statsOutput[thisRow + 2,] <- c(levels(pepisode$FrequencyBand)[thisFreq], "Pre1 vs Tele", pre1VsTele$p.value)
  statsOutput[thisRow + 3,] <- c(levels(pepisode$FrequencyBand)[thisFreq], "Tele vs Post1", teleVsPost1$p.value)
  statsOutput[thisRow + 4,] <- c(levels(pepisode$FrequencyBand)[thisFreq], "Post1 vs Post2", post1VsPost2$p.value)
  statsOutput[thisRow + 5,] <- c(levels(pepisode$FrequencyBand)[thisFreq], "Post2 vs Post3", post2VsPost3$p.value)
  
  thisRow <- thisRow + 6
}



