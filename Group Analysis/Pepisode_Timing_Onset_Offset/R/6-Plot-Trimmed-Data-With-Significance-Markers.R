# Script Name:  6-Plot-Trimmed-Data-With-Significance-Markers.R
# Author:       Lindsay Vass
# Date:         12 June 2015
# Purpose:      This script will recreate the plots in 4b-Plot-Trimmed-Data.R 
#               but will add significance markers based on the analyses in
#               5-Compute-Statistics-Trimmed-Data.

library(ggplot2)
library(grid)
library(ggthemes)
library(scales)
library(dplyr)

load('Rda/allTrimmedData.Rda')
load('Rda/allStats.Rda')

dir.create("Figures/TrimmedWithSigMarkers")

# Functions ---------------------------------------------------------------

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
makeGroupHistogram <- function(inputData, colour, fill, binBreaks, title, timePointMarkers, chiSquareMarkers, sigMarkers) {
  
  # get expected frequency based on number of observations
#   expFreq <- inputData %>%
#     group_by(ObservationType) %>%
#     summarise(Freq = n() / length(binBreaks))
  
  plotData <- inputData %>%
    ggplot(aes(x = Time)) +
    stat_bin(data = inputData,
             colour = colour,
             fill = fill,
             right = TRUE,
             breaks = binBreaks) +
    facet_grid(ObservationType ~ .) +
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
  
  sigMarkers <- inner_join(sigMarkers, facetInfo)
  sigMarkers <- inner_join(sigMarkers, histInfo)
  
  # positioning for statistics label
  yMaxInfo <- data.frame(PANEL = NA, y = NA)
  for (thisRow in 1:nrow(facetInfo)) {
    yMaxInfo[nrow(yMaxInfo) + 1, ] <- c(PANEL = thisRow,
                                        y = plotInfo$panel$ranges[[thisRow]]$y.range[2] +
                                          0.1 * plotInfo$panel$ranges[[thisRow]]$y.range[2])
  }
  yMaxInfo <- yMaxInfo %>% filter(is.na(PANEL) == FALSE)
  yMaxInfo$PANEL <- factor(yMaxInfo$PANEL)
  
  xInfo <- data.frame(PANEL = NA, x = NA)
  for (thisRow in 1:nrow(facetInfo)) {
    xInfo[nrow(xInfo) + 1, ] <- c(PANEL = thisRow,
                                  x = plotInfo$panel$ranges[[thisRow]]$x.range[1])
  }
  xInfo <- xInfo %>% filter(is.na(PANEL) == FALSE)
  xInfo$PANEL <- factor(xInfo$PANEL)
  
  
  facetInfo <- inner_join(facetInfo, yMaxInfo)
  facetInfo <- inner_join(facetInfo, xInfo)
  chiSquareMarkers <- inner_join(chiSquareMarkers, facetInfo)
  
  if (nrow(sigMarkers) == 0) {
    labeledPlotData <- plotData + 
      geom_text(data = chiSquareMarkers, aes(x = x, y = y, label = PlotLabel, hjust = 0))
  } else {
    labeledPlotData <- plotData + 
      geom_text(data = chiSquareMarkers, aes(x = x, y = y, label = PlotLabel, hjust = 0)) +
      geom_text(data = sigMarkers, size = 14, aes(x = x, y = y, label = "*", vjust = -0.1))
  }
  
}


# Make onset and offset histogram for each frequency/time type ------------

# These will be used to add vertical lines to the plots to indicate the times of
# teleporter entry and exit
timepointMarkers <- data.frame(x = c(0, 1830, 0, 2830), 
                               TrialTimeType = c("NT", "NT", "FT", "FT"))

# write chi-squared stats to label the plots
chiSquareMarkers <- chiSquareOutput %>%
  mutate(PlotLabel = paste0("Uncorrected P = ", signif(chiSquareOutput$RawPValue, digits = 2), 
                            "\nBonferroni Corrected (", nrow(chiSquareOutput), ") P = ", signif(chiSquareOutput$CorrectedPValue, digits = 2))) %>%
  dplyr::rename(ObservationType = Event) %>%
  select(-(ChiSquared:CorrectedPValue)) 

# convert to factors
frequencyBandOrder <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
timeTypeOrder      <- c("NT", "FT")
timepointMarkers$TrialTimeType <- factor(timepointMarkers$TrialTimeType, levels = timeTypeOrder)
chiSquareMarkers$FrequencyBand       <- factor(chiSquareMarkers$FrequencyBand, levels = frequencyBandOrder)
chiSquareMarkers$TrialTimeType       <- factor(chiSquareMarkers$TrialTimeType, levels = timeTypeOrder)
chiSquarePostHocOutput$FrequencyBand <- factor(chiSquarePostHocOutput$FrequencyBand, levels = frequencyBandOrder)
chiSquarePostHocOutput$TrialTimeType <- factor(chiSquarePostHocOutput$TrialTimeType, levels = timeTypeOrder)
chiSquarePostHocOutput$CorrectedPValue <- as.numeric(chiSquarePostHocOutput$CorrectedPValue)
chiSquarePostHocOutput <- dplyr::rename(chiSquarePostHocOutput, ObservationType = Event)

# set histogram parameters
colour     <- "steelblue4"
fill       <- "lightskyblue"
ntBinWidth <- 305 # time in ms, divides into 1830 evenly
ftBinWidth <- 283 # time in ms, divides into 2830 evenly
ntBinBreaks <- seq(-1830, 3660, ntBinWidth)
ftBinBreaks <- seq(-2830, 5660, ftBinWidth)

for (thisFreqBand in 1:nlevels(trimmedOnOffData$FrequencyBand)) {
  
  for (thisTimeType in 1:nlevels(trimmedOnOffData$TrialTimeType)) {
    
    # input data for the histogram
    thisData <- trimmedOnOffData %>%
      filter(FrequencyBand == levels(FrequencyBand)[thisFreqBand] &
               TrialTimeType == levels(TrialTimeType)[thisTimeType])
    
    # vertical lines to mark teleporter entry and exit
    theTimeMarkers <- timepointMarkers %>%
      filter(TrialTimeType == levels(TrialTimeType)[thisTimeType])
    
    # chi-squared statistics
    theChiSquareMarkers <- chiSquareMarkers %>%
      filter(FrequencyBand == levels(FrequencyBand)[thisFreqBand] &
               TrialTimeType == levels(TrialTimeType)[thisTimeType])
    
    # post-hoc binomial test results
    theSigMarkers <- chiSquarePostHocOutput %>%
      filter(FrequencyBand == levels(FrequencyBand)[thisFreqBand] &
               TrialTimeType == levels(TrialTimeType)[thisTimeType] &
               CorrectedPValue < 0.05)
    
    # size of bins
    if (levels(thisData$TrialTimeType)[thisTimeType] == "NT"){
      binBreaks <- ntBinBreaks
    } else {
      binBreaks <- ftBinBreaks
    }
                       
    # make histograms
    onsetPlot <- makeGroupHistogram(inputData = thisData,
                                    colour = colour,
                                    fill = fill,
                                    binBreaks = binBreaks,
                                    title = paste(levels(thisData$FrequencyBand)[thisFreqBand],
                                                  levels(thisData$TrialTimeType)[thisTimeType],
                                                  "Pepisode"),
                                    timePointMarkers = theTimeMarkers,
                                    chiSquareMarkers = theChiSquareMarkers,
                                    sigMarkers = theSigMarkers)
    ggsave(filename = paste0('Figures/TrimmedWithSigMarkers/', 
                             levels(thisData$FrequencyBand)[thisFreqBand],
                             '_',
                             levels(thisData$TrialTimeType)[thisTimeType],
                             '_Pepisode_Onset_Offset_Histograms.png'))
  }
}
