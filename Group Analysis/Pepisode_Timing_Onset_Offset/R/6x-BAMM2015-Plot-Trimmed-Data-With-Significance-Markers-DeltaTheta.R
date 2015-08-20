# Script Name:  6x-BAMM2015-Plot-Trimmed-Data-With-Significance-Markers.R
# Author:       Lindsay Vass
# Date:         18 August 2015
# Purpose:      This script will plot offset histograms for delta and theta for
#               presentation at BAMM 2015.

library(ggplot2)
library(grid)
library(ggthemes)
library(scales)
library(dplyr)

load('Rda/allTrimmedData.Rda')
load('Rda/allStats.Rda')

# Functions ---------------------------------------------------------------

# make strip text labels for scatterplots
scatterLabeller <- function(variable, value) {
  return(thisLabel$PlotLabel[value])
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

# the following two functions allow you to scale the y axis without specifying
# the limits
scale_dimension.custom_expand <- function(scale, expand = ggplot2:::scale_expand(scale)) {
  expand_range(ggplot2:::scale_limits(scale), expand[[1]], expand[[2]])
}
scale_y_continuous <- function(...) {
  s <- ggplot2::scale_y_continuous(...)
  class(s) <- c('custom_expand', class(s))
  s
}

# make histogram of group data
makeGroupHistogram <- function(inputData, colour, fill, binBreaks, title, timePointMarkers, sigMarkers) {
    
  plotData <- inputData %>%
    ggplot(aes(x = Time)) +
    stat_bin(data = inputData,
             colour = colour,
             fill = fill,
             right = TRUE,
             breaks = binBreaks) +
    geom_vline(aes(xintercept = x, size = 1), timePointMarkers) +
#     geom_hline(aes(yintercept = Freq), expFreq) +
    theme_few() +
    theme(text = element_text(size = 24),
          panel.margin = unit(1, "lines")) +
    labs(y = "# of Episodes", title = title) 
  
  # get info about the ylimits of each facet so we can add labels that appear
  # just below the max value and significance markers just above the bars
  plotInfo <- ggplot_build(plotData)
  facetInfo <- plotInfo$panel$layout
  histInfo <- plotInfo$data[[1]] %>%
    group_by(PANEL) 
  numBins <- histInfo %>% tally
  histInfo <- histInfo %>% 
    mutate(Bin = as.character(seq(1, numBins$n[1])))
  
  if (nrow(sigMarkers) != 0){
    sigMarkers <- inner_join(sigMarkers, facetInfo)
    sigMarkers <- inner_join(sigMarkers, histInfo)
  }
  
  
  # positioning for statistics label
#   yMaxInfo <- data.frame(PANEL = NA, y = NA)
#   for (thisRow in 1:nrow(facetInfo)) {
#     yMaxInfo[nrow(yMaxInfo) + 1, ] <- c(PANEL = thisRow,
#                                         y = plotInfo$panel$ranges[[thisRow]]$y.range[2] +
#                                           0.1 * plotInfo$panel$ranges[[thisRow]]$y.range[2])
#   }
#   yMaxInfo <- yMaxInfo %>% filter(is.na(PANEL) == FALSE)
#   yMaxInfo$PANEL <- factor(yMaxInfo$PANEL)
#   
#   xInfo <- data.frame(PANEL = NA, x = NA)
#   for (thisRow in 1:nrow(facetInfo)) {
#     xInfo[nrow(xInfo) + 1, ] <- c(PANEL = thisRow,
#                                   x = plotInfo$panel$ranges[[thisRow]]$x.range[1])
#   }
#   xInfo <- xInfo %>% filter(is.na(PANEL) == FALSE)
#   xInfo$PANEL <- factor(xInfo$PANEL)
#   
#   
#   facetInfo <- inner_join(facetInfo, yMaxInfo)
#   facetInfo <- inner_join(facetInfo, xInfo)
  
  if (nrow(sigMarkers) == 0) {
    labeledPlotData <- plotData 
  } else {
    labeledPlotData <- plotData + 
      geom_text(data = sigMarkers, size = 14, aes(x = x, y = y, label = "*", vjust = -0.1))
  }
  
}


# Make offset histogram for each frequency/time type ------------

# These will be used to add vertical lines to the plots to indicate the times of
# teleporter entry and exit
timepointMarkers <- data.frame(x = c(0, 1830, 0, 2830), 
                               TrialTimeType = c("NT", "NT", "FT", "FT"))

# convert to factors
frequencyBandOrder <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
timeTypeOrder      <- c("NT", "FT")
timepointMarkers$TrialTimeType <- factor(timepointMarkers$TrialTimeType, levels = timeTypeOrder)
chiSquarePostHocOutput$FrequencyBand <- factor(chiSquarePostHocOutput$FrequencyBand, levels = frequencyBandOrder)
chiSquarePostHocOutput$TrialTimeType <- factor(chiSquarePostHocOutput$TrialTimeType, levels = timeTypeOrder)
chiSquarePostHocOutput$CorrectedPValue <- as.numeric(chiSquarePostHocOutput$CorrectedPValue)
chiSquarePostHocOutput <- dplyr::rename(chiSquarePostHocOutput, ObservationType = Event)

# set histogram parameters
colour     <- "slategrey"
fill       <- "slategray2"
ntBinWidth <- 305 # time in ms, divides into 1830 evenly
ftBinWidth <- 283 # time in ms, divides into 2830 evenly
ntBinBreaks <- seq(-1830, 3660, ntBinWidth)
ftBinBreaks <- seq(-2830, 5660, ftBinWidth)

plotData <- trimmedOnOffData %>%
  filter(ObservationType == "Offset",
         FrequencyBand == "Delta" | FrequencyBand == "Theta") %>%
  group_by(FrequencyBand, TrialTimeType)


  for (thisTimeType in 1:nlevels(trimmedOnOffData$TrialTimeType)) {
    for (thisElectrode in 1:nlevels(trimmedOnOffData$ElectrodeID)) {
      
    
    # input data for the histogram
    thisData <- plotData %>%
      filter(TrialTimeType == levels(TrialTimeType)[thisTimeType],
             ElectrodeID == levels(ElectrodeID)[thisElectrode])
    
    # vertical lines to mark teleporter entry and exit
    theTimeMarkers <- timepointMarkers %>%
      filter(TrialTimeType == levels(TrialTimeType)[thisTimeType])
    
    # post-hoc binomial test results
    theSigMarkers <- chiSquarePostHocOutput %>%
      filter(TrialTimeType == levels(TrialTimeType)[thisTimeType],
             FrequencyBand == "Delta" | FrequencyBand == "Theta",
             ObservationType == "Offset",
             CorrectedPValue < 0.05)
    
    # size of histogram bins
    if (levels(thisData$TrialTimeType)[thisTimeType] == "NT"){
      binBreaks <- ntBinBreaks
    } else {
      binBreaks <- ftBinBreaks
    }
    
    # make histogram
    p <- thisData %>%
      ggplot(aes(x = Time)) +
      stat_bin(data = thisData,
               colour = colour,
               fill = fill,
               right = TRUE,
               breaks = binBreaks) +
      facet_grid(FrequencyBand ~ .) +
      geom_vline(aes(xintercept = x, size = 1), theTimeMarkers) +
      theme_stata() +
      theme(text = element_text(size = 24),
            axis.title.y = element_text(vjust = 1),
            panel.margin = unit(1, "lines"),
            plot.background = element_rect(fill = "white"),
            strip.background = element_rect(fill = fill)) +
      labs(y = "# of Episodes") 
    p
    
    
    
    
    
    ggsave(filename = paste0('Figures/BAMM2015/', 
                             levels(thisData$TrialTimeType)[thisTimeType],
                             '_Pepisode_Offset_Histograms_.png'),
           width = 20, height = 10)
    
  }
  
  }



for (thisFreqBand in 1:nlevels(trimmedOnOffData$FrequencyBand)) {
  
  for (thisTimeType in 1:nlevels(trimmedOnOffData$TrialTimeType)) {
    

                       
    # make histograms
    p <- makeGroupHistogram(inputData = thisData,
                            colour = colour,
                            fill = fill,
                            binBreaks = binBreaks,
                            title = paste(levels(thisData$FrequencyBand)[thisFreqBand],
                                          levels(thisData$TrialTimeType)[thisTimeType],
                                          "Pepisode"),
                            timePointMarkers = theTimeMarkers,
                            sigMarkers = theSigMarkers)
    
    
    
    
    ggsave(filename = paste0('Figures/TrimmedWithSigMarkers/', 
                             levels(thisData$FrequencyBand)[thisFreqBand],
                             '_',
                             levels(thisData$TrialTimeType)[thisTimeType],
                             '_Pepisode_Onset_Offset_Histograms.png'),
           width = 20, height = 10)
  }
}
