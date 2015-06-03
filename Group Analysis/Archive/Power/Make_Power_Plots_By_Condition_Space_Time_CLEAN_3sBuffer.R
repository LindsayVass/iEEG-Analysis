library(dplyr)
library(ggplot2)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power/')


## Load in the data for each session for each subject
expDir <- ('/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/')
subjectIDs <- c("UCDMC13","UCDMC14","UCDMC15")
teleporters <- c("TeleporterA","TeleporterB")

allPowerData <- data.frame()

for (thisSubject in 1:length(subjectIDs)) {
  
  for (thisTeleporter in 1:length(teleporters)) {
    
    # UCDMC13 doesn't have a TeleporterB, so skip it
    if (thisSubject == 1 & thisTeleporter == 2) {
      next
    }
    
    csvFile <- paste0(expDir, subjectIDs[thisSubject], "/Mat Files/Power/Summary/", subjectIDs[thisSubject], "_", teleporters[thisTeleporter], "_power_summary_noSpikes_noWaves_3sBuffer.csv")
    tempData <- read.csv(csvFile, header = TRUE)
    
    allPowerData <- rbind(allPowerData, tempData)
    
  } # thisTeleporter
  
} # thisSubject

## Put our conditions in order
timePointOrder <- c('Pre3', 'Pre2', 'Pre1', 'Tele', 'Post1', 'Post2', 'Post3')
allPowerData$TimePoint <- factor(allPowerData$TimePoint, levels = timePointOrder)

trialTypeOrder <- c('NSNT','NSFT','FSNT','FSFT')
allPowerData$TrialType <- factor(allPowerData$TrialType, levels = trialTypeOrder)

spaceOrder <- c('NS','FS')
allPowerData$TrialSpaceType <- factor(allPowerData$TrialSpaceType, levels = spaceOrder)

timeOrder <- c('NT','FT')
allPowerData$TrialTimeType <- factor(allPowerData$TrialTimeType, levels = timeOrder)

## Make a new variable for electrode ID (subjectID + Teleporter + Electrode)
allPowerData <- allPowerData %>%
  mutate(ElectrodeID = paste(SubjectID, Teleporter, Electrode, sep = "_")) %>%
  mutate(ElectrodeID = factor(ElectrodeID))

## Cut up the frequencies into bands
frequencies    <- unique(allPowerData$Frequency)
freqBandBreaks <- c(0, 4, 8, 12, 30, 182)
freqBandNames  <- c("Delta","Theta","Alpha","Beta","Gamma")

allPowerData <- allPowerData %>%
  mutate(FrequencyBand = cut(Frequency, freqBandBreaks, labels = freqBandNames))

# Write the data out for stats later
today <- Sys.Date()
format(today, format = "%d-%B-%Y")
saveFile <- paste0("All_Subjects_Power_CLEAN_3sBuffer", today, ".Rda")
save(allPowerData, file = saveFile)

##############################################
## summarize and plot the data by TIMEPOINT ##
##############################################
meanByTimePoint <- allPowerData %>%
  filter(is.na(Power) == FALSE) %>% # UCDMC15 TeleporterB is NaN for NSFT
  select(TimePoint,Power:FrequencyBand) %>%
  group_by(ElectrodeID, FrequencyBand,TimePoint) 

# We want our error bars to ultimately reflect the within-subject variability rather than between-subject
# To do this, we'll use the methods in Cousineau (2005) Tutorials in Quantitative Methods for Psychology and
# Morey (2008) [Same Journal], which corrects the bias in the Cousineau method
# Essentially, we will calculate normalized values for each observation, in which we take the original observation,
# subtract the electrode mean, add the group mean, and correct for the # of conditions

# calculate the mean across conditions for each electrode for each frequency
elecMean <- meanByTimePoint %>%
  group_by(ElectrodeID, FrequencyBand) %>%
  summarise(EMean = mean(Power))

# calculate the mean across conditions across electrodes for each frequency
groupMean <- meanByTimePoint %>%
  group_by(FrequencyBand) %>%
  summarise(GMean = mean(Power))

# now we'll normalize the observations by taking into account the electrode and group means
# to do this, we'll need to join our tables together
normByTimePoint <- inner_join(meanByTimePoint, elecMean)
normByTimePoint <- inner_join(normByTimePoint, groupMean) %>%
  mutate(NormPower = Power - EMean + GMean) %>%
  filter(is.na(Power) == FALSE) %>%
  group_by(FrequencyBand, TimePoint) %>%
  summarise(NormVar = var(NormPower)) %>%
  mutate(NormVarUnbias = NormVar * (length(timePointOrder)/(length(timePointOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(nlevels(allPowerData$ElectrodeID)))) %>%
  select(FrequencyBand, TimePoint, NormSEM)

# now summarise the original data and join with the normalized SEM
meanByTimePoint <- meanByTimePoint %>%
  group_by(FrequencyBand, TimePoint) %>%
  summarise(Power = mean(Power)) %>%
  inner_join(normByTimePoint)

# plot the data
ggplot(meanByTimePoint, aes(x = TimePoint, y = Power, ymin = Power - NormSEM, ymax = Power + NormSEM)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_wrap(~ FrequencyBand, scales = "free") + 
  ggtitle("Z-scored Log(Power) by Time Point")

ggsave(paste0('Figures/Power_by_timepoint_normSEM_3sBuffer', today, '.png'))


# make a plot showing Power by timepoint for each electrode
meanByElectrode <- allPowerData %>%
  filter(is.na(Power) == FALSE) %>%
  select(TimePoint, Power:FrequencyBand) %>%
  group_by(ElectrodeID, TimePoint, FrequencyBand) 
thetaByElectrode <- meanByElectrode %>%
  filter(FrequencyBand == "Theta") %>%
  summarise(PowerMean = mean(Power), SEM = sd(Power)/sqrt(n()))
ggplot(thetaByElectrode, aes(x = TimePoint, y = PowerMean, ymin = PowerMean - SEM, ymax = PowerMean + SEM)) +
  geom_point(size = 4) +
  geom_pointrange() + 
  facet_wrap(~ElectrodeID, scales = "free") +
  ggtitle("Theta Power by Electrode")


####################################################
## summarize and plot the data by SPACE condition ##
####################################################
meanBySpace <- allPowerData %>%
  filter(is.na(Power) == FALSE) %>% # UCDMC15 TeleporterB is NaN for EDF1 NSFT
  select(TrialSpaceType,TimePoint,Power:FrequencyBand) %>%
  group_by(ElectrodeID, FrequencyBand,TrialSpaceType) 

# We want our error bars to ultimately reflect the within-subject variability rather than between-subject
# To do this, we'll use the methods in Cousineau (2005) Tutorials in Quantitative Methods for Psychology and
# Morey (2008) [Same Journal], which corrects the bias in the Cousineau method
# Essentially, we will calculate normalized values for each observation, in which we take the original observation,
# subtract the electrode mean, add the group mean, and correct for the # of conditions

# calculate the mean across conditions for each electrode for each frequency
elecMean <- meanBySpace %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint, TrialSpaceType) %>%
  summarise(EMean = mean(Power))

# calculate the mean across conditions across electrodes for each frequency
groupMean <- meanBySpace %>%
  group_by(FrequencyBand, TimePoint, TrialSpaceType) %>%
  summarise(GMean = mean(Power))

# now we'll normalize the observations by taking into account the electrode and group means
# to do this, we'll need to join our tables together
normBySpace <- inner_join(meanBySpace, elecMean)
normBySpace <- inner_join(normBySpace, groupMean) %>%
  mutate(NormPower = Power - EMean + GMean) %>%
  group_by(FrequencyBand, TimePoint, TrialSpaceType) %>%
  summarise(NormVar = var(NormPower)) %>%
  mutate(NormVarUnbias = NormVar * (length(spaceOrder)/(length(spaceOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(nlevels(allPowerData$ElectrodeID)))) %>%
  select(FrequencyBand, TimePoint, TrialSpaceType, NormSEM)

# now summarise the original data and join with the normalized SEM
meanBySpace <- meanBySpace %>%
  group_by(FrequencyBand, TrialSpaceType, TimePoint) %>%
  summarise(Power = mean(Power)) %>%
  inner_join(normByTimePoint)

# plot the data
ggplot(meanBySpace, aes(x = TimePoint, y = Power, ymin = Power - NormSEM, ymax = Power + NormSEM, color = TrialSpaceType)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_wrap(~ FrequencyBand, scales = "free") + 
  ggtitle("Z-scored Log(Power) by Space Condition")

ggsave(paste0('Figures/Power_by_space_normSEM_3sBuffer', today, '.png'))

####################################################
## summarize and plot the data by Time condition ##
####################################################
meanByTime <- allPowerData %>%
  filter(is.na(Power) == FALSE) %>% # UCDMC15 TeleporterB is NaN for EDF1 NSFT
  select(TrialTimeType,TimePoint,Power:FrequencyBand) %>%
  group_by(ElectrodeID, FrequencyBand,TrialTimeType) 

# We want our error bars to ultimately reflect the within-subject variability rather than between-subject
# To do this, we'll use the methods in Cousineau (2005) Tutorials in Quantitative Methods for Psychology and
# Morey (2008) [Same Journal], which corrects the bias in the Cousineau method
# Essentially, we will calculate normalized values for each observation, in which we take the original observation,
# subtract the electrode mean, add the group mean, and correct for the # of conditions

# calculate the mean across conditions for each electrode for each frequency
elecMean <- meanByTime %>%
  group_by(ElectrodeID, FrequencyBand, TimePoint, TrialTimeType) %>%
  summarise(EMean = mean(Power))

# calculate the mean across conditions across electrodes for each frequency
groupMean <- meanByTime %>%
  group_by(FrequencyBand, TimePoint, TrialTimeType) %>%
  summarise(GMean = mean(Power))

# now we'll normalize the observations by taking into account the electrode and group means
# to do this, we'll need to join our tables together
normByTime <- inner_join(meanByTime, elecMean)
normByTime <- inner_join(normByTime, groupMean) %>%
  mutate(NormPower = Power - EMean + GMean) %>%
  group_by(FrequencyBand, TimePoint, TrialTimeType) %>%
  summarise(NormVar = var(NormPower)) %>%
  mutate(NormVarUnbias = NormVar * (length(timeOrder)/(length(timeOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(nlevels(allPowerData$ElectrodeID)))) %>%
  select(FrequencyBand, TimePoint, TrialTimeType, NormSEM)

# now summarise the original data and join with the normalized SEM
meanByTime <- meanByTime %>%
  group_by(FrequencyBand, TrialTimeType, TimePoint) %>%
  summarise(Power = mean(Power)) %>%
  inner_join(normByTimePoint)

# plot the data
ggplot(meanByTime, aes(x = TimePoint, y = Power, ymin = Power - NormSEM, ymax = Power + NormSEM, color = TrialTimeType)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_wrap(~ FrequencyBand, scales = "free") + 
  ggtitle("Z-scored Log(Power) by Time Condition")

ggsave(paste0('Figures/Power_by_time_normSEM_3sBuffer', today, '.png'))
