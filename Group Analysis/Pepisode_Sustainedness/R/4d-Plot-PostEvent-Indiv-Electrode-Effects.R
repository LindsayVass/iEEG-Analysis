# Date:     14 September 2015
# Purpose:  This script will plot the difference in post-teleporter-entry 
#           oscillatory event duration for navigation vs teleporter.

library(dplyr)
library(ggplot2)
library(ggthemes)

load('Rda/allAnalyzedData.Rda')
load('Rda/allCleanData.Rda')

# Functions ---------------------------------------------------------------

facetLabeller <- function(variable, value) {
  if (variable == "TimeType") {
    return(facetTimeLabels[value])
  } else{
    return(facetFreqLabels[value])
  }
}

meanPostEventOscDuration <- function(dataFrame, eventTime, trialTimeType) {
  output <- dataFrame %>%
    filter(Onset < eventTime & Offset > eventTime & TrialTimeType == trialTimeType) %>%
    mutate(PostEventDuration = Offset - eventTime) %>%
    group_by(ElectrodeID, FrequencyBand, RealTrialNumber) %>%
    summarise(MeanPostEventDuration = mean(PostEventDuration))
}


# Analyze data ------------------------------------------------------------

navNtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Navigation")
navFtPostEntryDur <- meanPostEventOscDuration(navSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Navigation")

teleNtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "NT") %>%
  mutate(TimeType = "NT",
         Condition = "Teleportation")
teleFtPostEntryDur <- meanPostEventOscDuration(teleSustain, 0, "FT") %>%
  mutate(TimeType = "FT",
         Condition = "Teleportation")
allPostEntryData <- rbind(navNtPostEntryDur, 
                          navFtPostEntryDur, 
                          teleNtPostEntryDur, 
                          teleFtPostEntryDur) %>%
  group_by(ElectrodeID, FrequencyBand, TimeType, Condition) %>%
  filter(n() > 5) %>%
  summarise(MeanDuration = mean(MeanPostEventDuration),
            SEM = sd(MeanPostEventDuration) / sqrt(n()))
navData <- allPostEntryData %>%
  filter(Condition == "Navigation")
teleData <- allPostEntryData %>%
  filter(Condition == "Teleportation")

validData <- inner_join(navData, teleData, by = c('ElectrodeID', 'FrequencyBand', 'TimeType')) %>%
  mutate(Difference = (MeanDuration.y - MeanDuration.x)) %>%
  select(ElectrodeID, TimeType, Difference, FrequencyBand) %>%
  inner_join(allPostEntryData)


# Make plot of individual electrode effects -------------------------------

validData$TimeType <- factor(validData$TimeType, c('NT', 'FT'))
facetTimeLabels <- c('Short Time', 'Long Time')
facetFreqLabels <- c("Delta-Theta", "Alpha", "Beta", "Gamma")
colFun <- colorRampPalette(c("red", "orange", "black", "deepskyblue", "dodgerblue4"))

p <- validData %>%
  ggplot(aes(x = Condition, 
             y = MeanDuration, 
             ymin = MeanDuration - SEM, 
             ymax = MeanDuration + SEM,
             group = ElectrodeID,
             colour = Difference)) +
  geom_point(size = 5) +
  geom_pointrange() +
  geom_line() +
  scale_color_gradientn(colours = colFun(5)) +
  theme_stata() +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 24),
        legend.text = element_blank(),
        legend.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1.5),
        strip.background = element_rect(colour = "black", size = 0.75),
        panel.border = element_rect(colour = "black", size = 0.75, fill = NA),
        panel.grid.major.y = element_line(colour = "dimgray", linetype = "longdash")) +
  ylab("Mean Duration (ms)") +
  facet_grid(TimeType ~ FrequencyBand, labeller = facetLabeller)
#dir.create('Figures/SingleElectrodeEventDur')
ggsave('Figures/SingleElectrodeEventDur/PostEntryDuration_Scatter.pdf', useDingbats = FALSE)