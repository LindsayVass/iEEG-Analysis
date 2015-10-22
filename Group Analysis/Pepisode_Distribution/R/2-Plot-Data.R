# Script Name:  2-Plot-Data.R
# Author:       Lindsay Vass
# Date:         15 October 2015
# Purpose:      This script will plot the distributions of power and pepisode
#               for each electrode

library(dplyr)
library(ggplot2)
library(ggthemes)
load('Rda/allRawData.Rda')

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  #   l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("^(.*)e", "e", l)
  # turn the 'e+' into plotmath format
  #   l <- gsub("e", "%*%10^", l)
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

pepisodeDf <- pepisodeData %>%
  mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = '_')) %>%
  group_by(ElectrodeID, Frequency) %>%
  select(-c(SubjectID, Teleporter, Electrode))
pepisodeDf$ElectrodeID <- factor(pepisodeDf$ElectrodeID)

powerDf <- powerData %>%
  mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = '_')) %>%
  group_by(ElectrodeID, Frequency) %>%
  select(-c(SubjectID, Teleporter, Electrode))
powerDf$ElectrodeID <- factor(powerDf$ElectrodeID)

freqLabels = c(2,4,8,16,32,64,128)

# Pepisode Bar Plot -------------------------------------------------------
for (i in 1:nlevels(pepisodeDf$ElectrodeID)) {
  thisData <- pepisodeDf %>%
    filter(ElectrodeID == levels(ElectrodeID)[i])
  p1 <- ggplot(thisData, aes(x = Frequency, y = Pepisode)) +
    geom_bar(stat = "identity") +
    scale_x_log10(breaks = freqLabels) +
    labs(x = 'Frequency (Hz)',
         y = expression('Mean P'['Episode'])) +
    theme_few() +
    theme(text = element_text(size = 24))
  ggsave(paste0('Figures/Pepisode_', levels(thisData$ElectrodeID)[i], '.pdf'))
}


# Power Spectral Density Plot ---------------------------------------------

for (i in 1:nlevels(powerDf$ElectrodeID)) {
  thisData <- powerDf %>%
    filter(ElectrodeID == levels(ElectrodeID)[i])
  p2 <- ggplot(thisData, aes(x = Frequency, y = Power)) +
    geom_line() +
    scale_x_log10(breaks = freqLabels) +
    scale_y_log10(labels = fancy_scientific) +
    labs(x = 'Frequency (Hz)',
         y = expression(paste('Power (', mu, 'V'^'2',' / Hz)'))) +
    theme_few() +
    theme(text = element_text(size = 24))
  ggsave(paste0('Figures/Power_', levels(thisData$ElectrodeID)[i], '.pdf'))
}

