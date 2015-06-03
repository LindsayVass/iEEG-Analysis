library(dplyr)
library(ggplot2)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/')


## Load in the data for each session for each subject
expDir <- ('/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/')
subjectIDs <- c("UCDMC13","UCDMC14","UCDMC15")
teleporters <- c("TeleporterA","TeleporterB")

allPepisodeData <- data.frame()

for (thisSubject in 1:length(subjectIDs)) {
  
  for (thisTeleporter in 1:length(teleporters)) {
    
    # UCDMC13 doesn't have a TeleporterB, so skip it
    if (thisSubject == 1 & thisTeleporter == 2) {
      next
    }
    
    csvFile <- paste0(expDir, subjectIDs[thisSubject], "/Mat Files/Pepisode/Summary/", subjectIDs[thisSubject], "_", teleporters[thisTeleporter], "_pepisode_summary_noSpikes_noWaves_TIMESHIFT.csv")
    tempData <- read.csv(csvFile, header = TRUE)
    
    allPepisodeData <- rbind(allPepisodeData, tempData)
    
  } # thisTeleporter
  
} # thisSubject

## Put our conditions in order
timePointOrder <- c('Pre','Tele','Post')
allPepisodeData$TimePoint <- factor(allPepisodeData$TimePoint, levels = timePointOrder)

trialTypeOrder <- c('NSNT','NSFT','FSNT','FSFT')
allPepisodeData$TrialType <- factor(allPepisodeData$TrialType, levels = trialTypeOrder)

spaceOrder <- c('NS','FS')
allPepisodeData$TrialSpaceType <- factor(allPepisodeData$TrialSpaceType, levels = spaceOrder)

timeOrder <- c('NT','FT')
allPepisodeData$TrialTimeType <- factor(allPepisodeData$TrialTimeType, levels = timeOrder)

## Make a new variable for electrode ID (subjectID + Teleporter + Electrode)
allPepisodeData <- allPepisodeData %>%
  mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = "_")) %>%
  mutate(ElectrodeID = factor(ElectrodeID))

## Cut up the frequencies into bands
frequencies    <- unique(allPepisodeData$Frequency)
freqBandBreaks <- c(0, 4, 8, 12, 30, 182)
freqBandNames  <- c("Delta","Theta","Alpha","Beta","Gamma")

allPepisodeData <- allPepisodeData %>%
  mutate(FrequencyBand = cut(Frequency, freqBandBreaks, labels = freqBandNames))

# Write the data out for stats later
write.csv(allPepisodeData, file = "All_Subjects_Pepisode_21May2015_CLEAN_TIMESHIFT.csv", row.names = FALSE)

##############################################
## summarize and plot the data by TIMEPOINT ##
##############################################
meanByTimePoint <- allPepisodeData %>%
  filter(is.na(Pepisode) == FALSE) %>% # UCDMC15 TeleporterB is NaN for NSFT
  select(TimePoint,Pepisode:FrequencyBand) %>%
  group_by(ElectrodeID, FrequencyBand,TimePoint) 

# We want our error bars to ultimately reflect the within-subject variability rather than between-subject
# To do this, we'll use the methods in Cousineau (2005) Tutorials in Quantitative Methods for Psychology and
# Morey (2008) [Same Journal], which corrects the bias in the Cousineau method
# Essentially, we will calculate normalized values for each observation, in which we take the original observation,
# subtract the electrode mean, add the group mean, and correct for the # of conditions

# calculate the mean across conditions for each electrode for each frequency
elecMean <- meanByTimePoint %>%
  group_by(ElectrodeID, FrequencyBand) %>%
  summarise(EMean = mean(Pepisode))

# calculate the mean across conditions across electrodes for each frequency
groupMean <- meanByTimePoint %>%
  group_by(FrequencyBand) %>%
  summarise(GMean = mean(Pepisode))

# now we'll normalize the observations by taking into account the electrode and group means
# to do this, we'll need to join our tables together
normByTimePoint <- inner_join(meanByTimePoint, elecMean)
normByTimePoint <- inner_join(normByTimePoint, groupMean) %>%
  mutate(NormPepisode = Pepisode - EMean + GMean) %>%
  filter(is.na(Pepisode) == FALSE) %>%
  group_by(FrequencyBand, TimePoint) %>%
  summarise(NormVar = var(NormPepisode)) %>%
  mutate(NormVarUnbias = NormVar * (length(timePointOrder)/(length(timePointOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(nlevels(allPepisodeData$ElectrodeID)))) %>%
  select(FrequencyBand, TimePoint, NormSEM)

# now summarise the original data and join with the normalized SEM
meanByTimePoint <- meanByTimePoint %>%
  group_by(FrequencyBand, TimePoint) %>%
  summarise(Pepisode = mean(Pepisode)) %>%
  inner_join(normByTimePoint)

# plot the data
ggplot(meanByTimePoint, aes(x = TimePoint, y = Pepisode, ymin = Pepisode - NormSEM, ymax = Pepisode + NormSEM)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_wrap(~ FrequencyBand, scales = "free") + 
  ggtitle("TIMESHIFT Pepisode by Time Point")

ggsave('Figures/TIMESHIFT_Pepisode_by_timepoint_normSEM_21May2015.png')


####################################################
## summarize and plot the data by SPACE condition ##
####################################################
meanBySpace <- allPepisodeData %>%
  filter(is.na(Pepisode) == FALSE) %>% # UCDMC15 TeleporterB is NaN for EDF1 NSFT
  select(TrialSpaceType,TimePoint,Pepisode:FrequencyBand) %>%
  group_by(ElectrodeID, FrequencyBand,TrialSpaceType) 

# We want our error bars to ultimately reflect the within-subject variability rather than between-subject
# To do this, we'll use the methods in Cousineau (2005) Tutorials in Quantitative Methods for Psychology and
# Morey (2008) [Same Journal], which corrects the bias in the Cousineau method
# Essentially, we will calculate normalized values for each observation, in which we take the original observation,
# subtract the electrode mean, add the group mean, and correct for the # of conditions

# calculate the mean across conditions for each electrode for each frequency
elecMean <- meanBySpace %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint, TrialSpaceType) %>%
  summarise(EMean = mean(Pepisode))

# calculate the mean across conditions across electrodes for each frequency
groupMean <- meanBySpace %>%
  group_by(FrequencyBand, TimePoint, TrialSpaceType) %>%
  summarise(GMean = mean(Pepisode))

# now we'll normalize the observations by taking into account the electrode and group means
# to do this, we'll need to join our tables together
normBySpace <- inner_join(meanBySpace, elecMean)
normBySpace <- inner_join(normBySpace, groupMean) %>%
  mutate(NormPepisode = Pepisode - EMean + GMean) %>%
  group_by(FrequencyBand, TimePoint, TrialSpaceType) %>%
  summarise(NormVar = var(NormPepisode)) %>%
  mutate(NormVarUnbias = NormVar * (length(spaceOrder)/(length(spaceOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(nlevels(allPepisodeData$ElectrodeID)))) %>%
  select(FrequencyBand, TimePoint, TrialSpaceType, NormSEM)

# now summarise the original data and join with the normalized SEM
meanBySpace <- meanBySpace %>%
  group_by(FrequencyBand, TrialSpaceType, TimePoint) %>%
  summarise(Pepisode = mean(Pepisode)) %>%
  inner_join(normByTimePoint)

# plot the data
ggplot(meanBySpace, aes(x = TimePoint, y = Pepisode, ymin = Pepisode - NormSEM, ymax = Pepisode + NormSEM, color = TrialSpaceType)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_wrap(~ FrequencyBand, scales = "free") + 
  ggtitle("TIMESHIFT Pepisode by Space Condition")

ggsave('Figures/TIMESHIFT_Pepisode_by_space_normSEM_21May2015.png')


####################################################
## summarize and plot the data by Time condition ##
####################################################
meanByTime <- allPepisodeData %>%
  filter(is.na(Pepisode) == FALSE) %>% # UCDMC15 TeleporterB is NaN for EDF1 NSFT
  select(TrialTimeType,TimePoint,Pepisode:FrequencyBand) %>%
  group_by(ElectrodeID, FrequencyBand,TrialTimeType) 

# We want our error bars to ultimately reflect the within-subject variability rather than between-subject
# To do this, we'll use the methods in Cousineau (2005) Tutorials in Quantitative Methods for Psychology and
# Morey (2008) [Same Journal], which corrects the bias in the Cousineau method
# Essentially, we will calculate normalized values for each observation, in which we take the original observation,
# subtract the electrode mean, add the group mean, and correct for the # of conditions

# calculate the mean across conditions for each electrode for each frequency
elecMean <- meanByTime %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint, TrialTimeType) %>%
  summarise(EMean = mean(Pepisode))

# calculate the mean across conditions across electrodes for each frequency
groupMean <- meanByTime %>%
  group_by(FrequencyBand, TimePoint, TrialTimeType) %>%
  summarise(GMean = mean(Pepisode))

# now we'll normalize the observations by taking into account the electrode and group means
# to do this, we'll need to join our tables together
normByTime <- inner_join(meanByTime, elecMean)
normByTime <- inner_join(normByTime, groupMean) %>%
  mutate(NormPepisode = Pepisode - EMean + GMean) %>%
  group_by(FrequencyBand, TimePoint, TrialTimeType) %>%
  summarise(NormVar = var(NormPepisode)) %>%
  mutate(NormVarUnbias = NormVar * (length(timeOrder)/(length(timeOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(nlevels(allPepisodeData$ElectrodeID)))) %>%
  select(FrequencyBand, TimePoint, TrialTimeType, NormSEM)

# now summarise the original data and join with the normalized SEM
meanByTime <- meanByTime %>%
  group_by(FrequencyBand, TrialTimeType, TimePoint) %>%
  summarise(Pepisode = mean(Pepisode)) %>%
  inner_join(normByTimePoint)

# plot the data
ggplot(meanByTime, aes(x = TimePoint, y = Pepisode, ymin = Pepisode - NormSEM, ymax = Pepisode + NormSEM, color = TrialTimeType)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_wrap(~ FrequencyBand, scales = "free") + 
  ggtitle("TIMESHIFT Pepisode by Time Condition")

ggsave('Figures/TIMESHIFT_Pepisode_by_Time_normSEM_21May2015.png')
