# Script Name:  4-Plot-Data.R
# Author:       Lindsay Vass
# Date:         18 June 2015
# Purpose:      This script will plot the data in
#               3b-Analyze-Data-DeltaThetaCombined.R. It will produce scatter 
#               plots for each time type showing mean pepisode for each time 
#               point (Pre/Tele/Post) and space type (NS/FS).

library(dplyr)
library(ggplot2)
library(grid)
library(ggthemes)
library(lazyeval)
library(reshape2)

load('Rda/allAnalyzedData_DeltaThetaCombined.Rda')
load('Rda/allCleanData.Rda')
dir.create('Figures')

# Functions ---------------------------------------------------------------

# return data from a specific frequency band and trial time type
filterRMData <- function(episodeData, timeType) {
  filteredData <- episodeData %>%
    filter(TrialTimeType == timeType) %>%
    ungroup() %>%
    group_by(ElectrodeID, TimeBin)
}

calcGroupMeanElecMean <- function(dataFrame, elecVar = "ElectrodeID", facetVar2 = "TrialTimeType", facetVar3 = "TrialSpaceType", conditionVar = "TimePoint", dataVar) {
  # We want our error bars to ultimately reflect the within-subject variability
  # rather than between-subject To do this, we'll use the methods in Cousineau
  # (2005) Tutorials in Quantitative Methods for Psychology and Morey (2008)
  # [Same Journal], which corrects the bias in the Cousineau method Essentially,
  # we will calculate normalized values for each observation, in which we take
  # the original observation, subtract the electrode mean, add the group mean,
  # and correct for the # of conditions
  
  # calculate the mean across conditions for each electrode grouped by facetVar
  elecMean <- dataFrame %>%
    group_by_(elecVar, facetVar2, facetVar3) %>%
    summarise_(EMean = interp(~mean(dataVar), dataVar = as.name(dataVar)))
  
  # calculate the mean across conditions across electrodes grouped by groupVar
  groupMean <- dataFrame %>%
    group_by_(facetVar2, facetVar3) %>%
    summarise_(GMean = interp(~mean(dataVar), dataVar = as.name(dataVar)))
  
  # now we'll normalize the observations by taking into account the electrode
  # and group means to do this, we'll need to join our tables together
  normData <- inner_join(dataFrame, elecMean)
  normData <- inner_join(normData, groupMean) %>%
    mutate_(NormData = interp(~(dataVar - EMean + GMean), dataVar = as.name(dataVar))) %>%
    filter_(interp(~(is.na(dataVar) == FALSE), dataVar = as.name(dataVar))) %>%
    group_by_(facetVar2, facetVar3, conditionVar) %>%
    summarise_(NormVar = interp(~var(NormData), NormData = as.name("NormData"))) %>%
    mutate_(NormVarUnbias = interp(~(NormVar * (nlevels(dataFrame$conditionVar) / (nlevels(dataFrame$conditionVar) - 1))), 
                                   conditionVar = as.name(conditionVar),
                                   NormVar = as.name("NormVar")),
            NormSEM = interp(~(sqrt(NormVarUnbias) / sqrt(nlevels(dataFrame$elecVar))), 
                             elecVar = as.name(elecVar),
                             NormVarUnbias = as.name("NormVarUnbias"))) 
  
  # now summarise the original data and join with the normalized SEM
  dataFrame <- dataFrame %>%
    group_by_(facetVar2, facetVar3, conditionVar) %>%
    summarise_(Value = interp(~mean(dataVar), dataVar = as.name(dataVar))) %>%
    inner_join(normData)
  
  return(dataFrame)
}

# make strip text labels for scatterplots
scatterLabeller <- function(variable, value) {
  return(thisLabel$PlotLabel[value])
}

# run posthoc test
runPosthoc <- function(inputData, timeType) {
  
  posthocData <- inputData %>%
    filter(TrialTimeType == timeType) %>%
    dcast(ElectrodeID ~ TrialSpaceType, value.var = "Pepisode")
  result <- wilcox.test(posthocData$NS, posthocData$FS, paired = TRUE)
}

# Prep data for scatterplots ----------------------------------------------

cleanData <- cleanData %>%
  filter(FrequencyBand == "Delta" | FrequencyBand == "Theta")

# get repeated-measures error bars
plotData <- calcGroupMeanElecMean(cleanData, 
                                  elecVar = "ElectrodeID", 
                                  facetVar2 = "TrialTimeType", 
                                  facetVar3 = "TrialSpaceType",
                                  conditionVar = "TimeBin", 
                                  dataVar = "Pepisode")

# make ANOVA text to print on plots
anovaResults <- anovaResults %>%
  group_by(TimeType) %>%
  mutate(ContrastLabel = paste0(Contrast,
                                ifelse(Contrast == "spaceType:timePoint", '\tF = ', '\t\tF = '),
                                sprintf("%.2f", Fstat),
                                '\tPerm. Corrected P = ',
                                sprintf("%.3f", Pcorr))) %>%
  dcast(TimeType ~ Contrast, value.var = "ContrastLabel") %>%
  rename(spaceTypeXtimePoint = `spaceType:timePoint`) %>% #paste0 can't interpret the original column name in the next step
  mutate(PlotLabel = paste0(TimeType,
                            '\n',
                            spaceType,
                            '\n',
                            timePoint,
                            '\n',
                            spaceTypeXtimePoint)) %>%
  select(TimeType, PlotLabel)
anovaResults$TimeType <- factor(anovaResults$TimeType, levels = levels(cleanData$TrialTimeType))


# plot the data 
thisData <- plotData 
thisLabel <- anovaResults %>%
  arrange(TimeType)

ggplot(thisData, aes(x = TimeBin, y = Value, ymin = Value - NormSEM, ymax = Value + NormSEM, color = TrialSpaceType)) +
  geom_point(size = 4) +
  geom_pointrange() +
  facet_grid(~TrialTimeType, labeller = scatterLabeller) +
  ggtitle("Delta/Theta Pepisode by Time Point") +
  theme(text = element_text(size = 18),
        strip.text = element_text(hjust = 0)) +
  ylab("Pepisode")
ggsave(filename = 'Figures/DeltaTheta_Pepisode_by_Time_Space_and_TimePoint.png',
       width = 16, height = 8)   

# posthoc tests on space x timepoint interaction
posthocData <- cleanData %>%
  filter(TimeBin == "Tele") %>%
  ungroup() %>%
  group_by(ElectrodeID, TrialSpaceType, TrialTimeType) %>%
  summarise(Pepisode = mean(Pepisode))
ntPosthoc <- runPosthoc(posthocData, "NT")
ftPosthoc <- runPosthoc(posthocData, "FT")

posthocResults <- data.frame(FrequencyBand = "Delta/Theta", 
                             TimePoint = "Tele", 
                             TimeType = c("NT", "FT"),
                             VStat = c(ntPosthoc$statistic, ftPosthoc$statistic),
                             PValue = c(ntPosthoc$p.value, ftPosthoc$p.value))

# save
save(file = 'Rda/allAnalyzedData_DeltaThetaCombined.Rda', list = c('anovaResults', 'posthocResults'))
  