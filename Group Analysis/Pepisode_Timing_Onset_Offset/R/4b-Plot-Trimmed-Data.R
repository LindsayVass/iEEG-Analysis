# Script Name:  4b-Plot-Trimmed-Data.R
# Author:       Lindsay Vass
# Date:         9 June 2015
# Purpose:      This script will plot the data produced by 3-Analyze-Data.R. It
#               will produce line plots that show when pepisodes are present
#               over the course of the trial as well as histograms that show
#               when episodes are most likely to begin and end. It differs from
#               4-Plot-Data.R in that the data are trimmed to remove edge effects.

library(ggplot2)
library(dplyr)
library(ggthemes)
library(grid)



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
    facet_grid(FrequencyBand ~ TrialTimeType, scales = "free") +
    geom_vline(aes(xintercept = x), timePointMarkers) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = title)
}

# make electrodewise histograms
makeElectrodewiseHistogram <- function(inputData, colour, fill, binWidth, title, timePointMarkers) {
  
  plotData <- inputData %>%
    ggplot(aes(x = Time)) +
    stat_bin(colour = colour,
             fill = fill,
             right = TRUE,
             binwidth = binWidth) +
    facet_wrap(~ElectrodeID, scales = "free") +
    geom_vline(xintercept = timePointMarkers[1]) +
    geom_vline(xintercept = timePointMarkers[2]) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = title)
  
}

# Trim the data -----------------------------------------------------------
load('Rda/allAnalyzedData.Rda')

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

# Make line plots ---------------------------------------------------------

# group data
timepointMarkers <- data.frame(x = c(0, 1830, 0, 2830), TrialTimeType = c("NT", "NT", "FT", "FT"))
lineData <- trimmedEpisodeData %>%
  group_by(FrequencyBand, TrialTimeType, Time) %>%
  summarise(GroupMean = mean(Mean), GroupSEM = sd(Mean) / sqrt(n())) 
linePlot <- lineData %>%
  ggplot(aes(x = Time, y = GroupMean, ymin = GroupMean - GroupSEM, ymax = GroupMean + GroupSEM)) + 
  geom_ribbon(color = "lightskyblue", fill = "lightskyblue") +
  geom_line(color = "steelblue4") +
  geom_vline(aes(xintercept = x, linetype = "dashed"), timepointMarkers) +
  facet_grid(FrequencyBand ~ TrialTimeType, scales = "free_x") +
  theme_few() +
  labs(y = "Mean Pepisode", title = "Pepisode Over Time") +
  theme(text = element_text(size = 24),
        axis.text = element_text(size = 18),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(vjust = 1),
        panel.margin = unit(1, "lines"))
dir.create("Figures/Trimmed")
ggsave('Figures/Trimmed/Pepisode_Over_Time.png', linePlot)

# individual electrode data
for (thisFrequency in 1:nlevels(episodeData$FrequencyBand)) {
  lineNtData <- trimmedEpisodeData %>%
    filter(FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency] &
             TrialTimeType == "NT")%>%
    ggplot(aes(x = Time, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM)) + 
    geom_ribbon(color = "lightskyblue", fill = "lightskyblue") +
    geom_line(color = "steelblue4") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 1830) +
    facet_wrap(~ElectrodeID) +
    theme_few() +
    labs(y = "Mean Pepisode", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Over Time (NT Trials)")) +
    theme(text = element_text(size = 24),
          axis.text = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 1),
          panel.margin = unit(1, "lines"))
  
    ggsave(filename = paste0('Figures/Trimmed/', levels(episodeData$FrequencyBand)[thisFrequency], '_NT_Electrodewise_Pepisode_Over_Time.png'),
           plot = lineNtData,
           width = 30,
           height = 15)

  lineFtData <- trimmedEpisodeData %>%
    filter(FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency] &
             TrialTimeType == "FT")%>%
    ggplot(aes(x = Time, y = Mean, ymin = Mean - SEM, ymax = Mean + SEM)) + 
    geom_ribbon(color = "lightskyblue", fill = "lightskyblue") +
    geom_line(color = "steelblue4") +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 2830) +
    facet_wrap(~ElectrodeID) +
    theme_few() +
    labs(y = "Mean Pepisode", title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], " Pepisode Over Time (FT Trials)")) +
    theme(text = element_text(size = 24),
          axis.text = element_text(size = 18),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 1),
          panel.margin = unit(1, "lines"))
  
  ggsave(filename = paste0('Figures/Trimmed/', levels(episodeData$FrequencyBand)[thisFrequency], '_FT_Electrodewise_Pepisode_Over_Time.png'),
         plot = lineFtData,
         width = 30,
         height = 15)
}

# Make onset histograms ------------------------------------------------------

# set histogram parameters
colour     <- "steelblue4"
fill       <- "lightskyblue"
ntBinWidth <- 305
ftBinWidth <- 283


# group data
onsetData <- trimmedOnOffData %>%
  ungroup() %>%
  filter(ObservationType == "Onset")
onsetPlot <- makeGroupHistogram(inputData = onsetData,
                                colour,
                                fill,
                                ntBinWidth,
                                ftBinWidth,
                                title = "Pepisode Onset Times",
                                timepointMarkers)
ggsave('Figures/Trimmed/Pepisode_Onset_Times_Histogram.png', onsetPlot)

# electrodewise data
for (thisFrequency in 1:nlevels(episodeData$FrequencyBand)) {
  onsetNtData <- trimmedOnOffData %>%
    ungroup() %>%
    filter(ObservationType == "Onset" &  
             TrialTimeType == "NT" &
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) 
  plotNtData <- makeElectrodewiseHistogram(inputData = onsetNtData,
                                          colour,
                                          fill,
                                          ntBinWidth,
                                          title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], 
                                                         " Pepisode Onset Times (NT Trials)"),
                                          c(0, 1830))
  ggsave(filename =  paste0('Figures/Trimmed/', levels(episodeData$FrequencyBand)[thisFrequency], '_NT_Electrodewise_Pepisode_Onset_Times_Histogram.png'),
         plot = plotNtData,
         width = 30,
         height = 15)
  
  onsetFtData <- trimmedOnOffData %>%
    ungroup() %>%
    filter(ObservationType == "Onset" &  
             TrialTimeType == "FT" &
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) 
  plotFtData <- makeElectrodewiseHistogram(inputData = onsetFtData,
                                           colour,
                                           fill,
                                           ftBinWidth,
                                           title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], 
                                                          " Pepisode Onset Times (FT Trials)"),
                                           c(0, 2830))
  ggsave(filename =  paste0('Figures/Trimmed/', levels(episodeData$FrequencyBand)[thisFrequency], '_FT_Electrodewise_Pepisode_Onset_Times_Histogram.png'),
         plot = plotFtData,
         width = 30,
         height = 15)
         
}
 

# Make offset histograms ------------------------------------------------------

# group data
offsetData <- trimmedOnOffData %>%
  ungroup() %>%
  filter(ObservationType == "Offset")
offsetPlot <- makeGroupHistogram(inputData = offsetData,
                                colour,
                                fill,
                                ntBinWidth,
                                ftBinWidth,
                                title = "Pepisode Offset Times",
                                timepointMarkers)
ggsave('Figures/Trimmed/Pepisode_Offset_Times_Histogram.png', offsetPlot)

# electrodewise data
for (thisFrequency in 1:nlevels(episodeData$FrequencyBand)) {
  offsetNtData <- trimmedOnOffData %>%
    ungroup() %>%
    filter(ObservationType == "Offset" &  
             TrialTimeType == "NT" &
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) 
  plotNtData <- makeElectrodewiseHistogram(inputData = offsetNtData,
                                           colour,
                                           fill,
                                           ntBinWidth,
                                           title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], 
                                                          " Pepisode Offset Times (NT Trials)"),
                                           c(0, 1830))
  ggsave(filename =  paste0('Figures/Trimmed/', levels(episodeData$FrequencyBand)[thisFrequency], '_NT_Electrodewise_Pepisode_Offset_Times_Histogram.png'),
         plot = plotNtData,
         width = 30,
         height = 15)
  
  offsetFtData <- trimmedOnOffData %>%
    ungroup() %>%
    filter(ObservationType == "Offset" &  
             TrialTimeType == "FT" &
             FrequencyBand == levels(episodeData$FrequencyBand)[thisFrequency]) 
  plotFtData <- makeElectrodewiseHistogram(inputData = offsetFtData,
                                           colour,
                                           fill,
                                           ftBinWidth,
                                           title = paste0(levels(episodeData$FrequencyBand)[thisFrequency], 
                                                          " Pepisode Offset Times (FT Trials)"),
                                           c(0, 2830))
  ggsave(filename =  paste0('Figures/Trimmed/', levels(episodeData$FrequencyBand)[thisFrequency], '_FT_Electrodewise_Pepisode_Offset_Times_Histogram.png'),
         plot = plotFtData,
         width = 30,
         height = 15)
  
}


# Plot mean pepisode by time bin ------------------------------------------

plotEpisodeData <- calcGroupMeanElecMean(binnedEpisodeData, 
                                         elecVar = "ElectrodeID", 
                                         facetVar1 = "FrequencyBand", 
                                         facetVar2 = "TrialTimeType", 
                                         conditionVar = "TimeBin", 
                                         dataVar = "Mean")


# plot the data
ggplot(plotEpisodeData, aes(x = TimeBin, y = Value, ymin = Value - NormSEM, ymax = Value + NormSEM)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_grid(FrequencyBand ~ TrialTimeType, scales = "free") + 
  ggtitle("Pepisode by Time Point")


# Save data ---------------------------------------------------------------

save(list = c("binnedEpisodeData", "trimmedOnOffData"), file = 'Rda/allTrimmedData.Rda')
