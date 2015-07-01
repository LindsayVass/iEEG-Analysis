# Script Name:  4c-Plot-Trimmed-Data-DeltaThetaEachFrequency.R
# Author:       Lindsay Vass
# Date:         25 June 2015
# Purpose:      This script will plot the data produced by 3-Analyze-Data.R. It
#               will produce histograms that show
#               when episodes are most likely to begin and end. It differs from
#               4-Plot-Data.R in that the data are trimmed to remove edge effects.
#               It differs from 4b-Plot-Data.R in that we will separately 
#               plot the data for each valid frequency within the delta and
#               theta bands.

library(ggplot2)
library(dplyr)
library(ggthemes)
library(grid)

load('Rda/allAnalyzedData.Rda')

# Functions ---------------------------------------------------------------

factors <- function(n) {
  if (length(n) > 1) {
    lapply(as.list(n), factors)
  } else {
    one.to.n <- seq_len(n)
    one.to.n[(n %% one.to.n) == 0]
  }
}

calcGroupMeanElecMean <- function(dataFrame, elecVar = "ElectrodeID", facetVar1 = "FrequencyBand", facetVar2 = "TrialTimeType", conditionVar = "TimePoint", dataVar) {
  # We want our error bars to ultimately reflect the within-subject variability
  # rather than between-subject To do this, we'll use the methods in Cousineau
  # (2005) Tutorials in Quantitative Methods for Psychology and Morey (2008)
  # [Same Journal], which corrects the bias in the Cousineau method Essentially,
  # we will calculate normalized values for each observation, in which we take
  # the original observation, subtract the electrode mean, add the group mean,
  # and correct for the # of conditions
  
  # calculate the mean across conditions for each electrode grouped by facetVar
  elecMean <- dataFrame %>%
    group_by_(elecVar, facetVar1, facetVar2) %>%
    summarise_(EMean = interp(~mean(dataVar), dataVar = as.name(dataVar)))
  
  # calculate the mean across conditions across electrodes grouped by groupVar
  groupMean <- dataFrame %>%
    group_by_(facetVar1, facetVar2) %>%
    summarise_(GMean = interp(~mean(dataVar), dataVar = as.name(dataVar)))
  
  # now we'll normalize the observations by taking into account the electrode
  # and group means to do this, we'll need to join our tables together
  normData <- inner_join(dataFrame, elecMean)
  normData <- inner_join(normData, groupMean) %>%
    mutate_(NormData = interp(~(dataVar - EMean + GMean), dataVar = as.name(dataVar))) %>%
    filter_(interp(~(is.na(dataVar) == FALSE), dataVar = as.name(dataVar))) %>%
    group_by_(facetVar1, facetVar2, conditionVar) %>%
    summarise_(NormVar = interp(~var(NormData), NormData = as.name("NormData"))) %>%
    mutate_(NormVarUnbias = interp(~(NormVar * (nlevels(dataFrame$conditionVar) / (nlevels(dataFrame$conditionVar) - 1))), 
                                   conditionVar = as.name(conditionVar),
                                   NormVar = as.name("NormVar")),
           NormSEM = interp(~(sqrt(NormVarUnbias) / sqrt(nlevels(dataFrame$elecVar))), 
                            elecVar = as.name(elecVar),
                            NormVarUnbias = as.name("NormVarUnbias"))) 
  
  # now summarise the original data and join with the normalized SEM
  dataFrame <- dataFrame %>%
    group_by_(facetVar1, facetVar2, conditionVar) %>%
    summarise_(Value = interp(~mean(dataVar), dataVar = as.name(dataVar))) %>%
    inner_join(normData)
  
  return(dataFrame)
}

# make histogram of group data
makeGroupHistogram <- function(inputData, colour, fill, ntBinWidth, ftBinWidth, title, timePointMarkers) {
  plotData <- inputData %>%
    ggplot(aes(x = Time)) +
    stat_bin(data = subset(inputData, TrialTimeType == "NT"),
             colour = colour,
             fill = fill,
             right = TRUE,
             binwidth = ntBinWidth) +
    stat_bin(data = subset(inputData, TrialTimeType == "FT"),
             colour = colour,
             fill = fill,
             right = TRUE,
             binwidth = ftBinWidth) +
    facet_grid(ObservationType ~ TrialTimeType, scales = "free_x") +
    geom_vline(aes(xintercept = x), timePointMarkers) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = title)
}

# Trim the data -----------------------------------------------------------

ntTrim <- c(-1830, 3660)
ftTrim <- c(-2830, 5660)

ntData <- episodeData %>%
  filter(TrialTimeType == "NT" &
           Time >= ntTrim[1] &
           Time <= ntTrim[2])
ftData <- episodeData %>%
  filter(TrialTimeType == "FT" &
           Time >= ftTrim[1] &
           Time <= ftTrim[2])
trimmedEpisodeData <- rbind(ntData, ftData)

# Divide time points into pre/tele/post bins
timeBinBreaksNT <- c(-1830, 0, 1830, 3660)
timeBinBreaksFT <- c(-2830, 0, 2830, 5660)
timeBinNames  <- c("Pre", "Tele", "Post")

binnedEpisodeData <- trimmedEpisodeData %>%
  mutate(TimeBin = if(TrialTimeType == "NT") {
    TimeBin = cut(Time, timeBinBreaksNT, labels = timeBinNames)
  } else {
    TimeBin = cut(Time, timeBinBreaksFT, labels = timeBinNames)
  })


ntData <- onOffData %>%
  ungroup() %>%
  filter(TrialTimeType == "NT" &
           Time >= ntTrim[1] &
           Time <= ntTrim[2])
ftData <- onOffData %>%
  ungroup() %>%
  filter(TrialTimeType == "FT" &
           Time >= ftTrim[1] &
           Time <= ftTrim[2])
trimmedOnOffData <- rbind(ntData, ftData)


# remove frequencies outside delta/theta
minFreq <- 1000 / (1830 / 3)
trimmedOnOffData <- trimmedOnOffData %>%
  filter(Frequency >= minFreq & Frequency <= 8)

trimmedOnOffData$Frequency <- factor(trimmedOnOffData$Frequency)
trimmedOnOffData$ObservationType <- factor(trimmedOnOffData$ObservationType)

# Make onset and offset histograms for each frequency ------------------------

# set histogram parameters
colour     <- "steelblue4"
fill       <- "lightskyblue"
ntBinWidth <- 305
ftBinWidth <- 283
timepointMarkers <- data.frame(x = c(0, 1830, 0, 2830), TrialTimeType = c("NT", "NT", "FT", "FT"))

for (thisFreq in 1:nlevels(trimmedOnOffData$Frequency)) {
  
  thisData <- trimmedOnOffData %>%
    ungroup() %>%
    filter(Frequency == levels(Frequency)[thisFreq])
  thisPlot <- makeGroupHistogram(thisData,
                                 colour,
                                 fill,
                                 ntBinWidth,
                                 ftBinWidth,
                                 title = paste(substr(levels(thisData$Frequency)[thisFreq], 1, 4), "Hz Pepisode Onset and Offset Times"),
                                 timepointMarkers)
  ggsave(paste0('Figures/Trimmed/', substr(levels(thisData$Frequency)[thisFreq], 1, 4), 'Hz_Pepisode_Onset_Offset_Times_Histogram.png'), thisPlot)
  
}


# Make onset and offset histograms collapsed across frequencies -----------

collapsedData <- trimmedOnOffData %>%
  ungroup()
collapsedPlot <- makeGroupHistogram(collapsedData,
                                    colour,
                                    fill,
                                    ntBinWidth,
                                    ftBinWidth,
                                    title = "Delta/Theta Pepisode Onset and Offset Times",
                                    timepointMarkers)
ggsave('Figures/Trimmed/DeltaTheta_Pepisode_Onset_Offset_Times_Histogram.png')

# Save data ---------------------------------------------------------------

save(list = c("trimmedOnOffData"), file = 'Rda/allTrimmedData_DeltaThetaByFrequency.Rda')
