library(dplyr)
library(ggplot2)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/')
filePath <- 'Tall_Pepisode_by_condition_30-Apr-2015.csv'
pepisode <- read.csv(filePath, header = FALSE)
names(pepisode) <- c('Electrode','Condition','FreqBand','TimePoint','Power')

# put our conditions in order
timePointOrder <- c('Pre','Tele','Post')
conditionOrder <- c('NSNT','NSFT','FSNT','FSFT')
freqOrder <- c('Delta','Theta','Alpha','Beta','Gamma')
pepisode$Condition <- factor(pepisode$Condition, levels = conditionOrder)
pepisode$TimePoint <- factor(pepisode$TimePoint, levels = timePointOrder)
pepisode$FreqBand <- factor(pepisode$FreqBand, levels = freqOrder)

##############################################
## summarize and plot the data by TIMEPOINT ##
##############################################
meanByTimePoint <- pepisode %>%
  select(-Condition) %>%
  group_by(Electrode, FreqBand,TimePoint) 

# We want our error bars to ultimately reflect the within-subject variability rather than between-subject
# To do this, we'll use the methods in Cousineau (2005) Tutorials in Quantitative Methods for Psychology and
# Morey (2008) [Same Journal], which corrects the bias in the Cousineau method
# Essentially, we will calculate normalized values for each observation, in which we take the original observation,
# subtract the electrode mean, add the group mean, and correct for the # of conditions

# calculate the mean across conditions for each electrode for each frequency
elecMean <- meanByTimePoint %>%
  group_by(Electrode, FreqBand) %>%
  summarise(EMean = mean(Power))

# calculate the mean across conditions across electrodes for each frequency
groupMean <- meanByTimePoint %>%
  group_by(FreqBand) %>%
  summarise(GMean = mean(Power))

# now we'll normalize the observations by taking into account the electrode and group means
# to do this, we'll need to join our tables together
normByTimePoint <- inner_join(meanByTimePoint, elecMean)
normByTimePoint <- inner_join(normByTimePoint, groupMean) %>%
  mutate(NormPepisode = Power - EMean + GMean) %>%
  group_by(FreqBand, TimePoint) %>%
  summarise(NormVar = var(NormPepisode)) %>%
  mutate(NormVarUnbias = NormVar * (length(timePointOrder)/(length(timePointOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(max(pepisode$Electrode)))) %>%
  select(FreqBand, TimePoint, NormSEM)

# now summarise the original data and join with the normalized SEM
meanByTimePoint <- meanByTimePoint %>%
  group_by(FreqBand, TimePoint) %>%
  summarise(Pepisode = mean(Power)) %>%
  inner_join(normByTimePoint)

# plot the data
ggplot(meanByTimePoint, aes(x = TimePoint, y = Pepisode, ymin = Pepisode - NormSEM, ymax = Pepisode + NormSEM)) +
  geom_point(size = 4) + 
  geom_pointrange() +
  facet_wrap(~ FreqBand, scales = "free") + 
  ggtitle("Pepisode by Time Point")

ggsave('Figures/Pepisode_by_timepoint_normSEM_11May2015.png')


####################################################
## summarize and plot the data by SPACE condition ##
####################################################
meanBySpace <- pepisode %>%
  mutate(Space = substr(as.character(Condition),1,2)) %>%
  select(-Condition) %>%
  group_by(Electrode, FreqBand,TimePoint)
spaceOrder <- c('NS','FS')
meanBySpace$Space <- factor(meanBySpace$Space, levels = spaceOrder)

# calculate mean across conditions for each electrode for each frequency
elecMean <- meanBySpace %>%
  group_by(Space, Electrode, FreqBand) %>%
  summarise(EMean = mean(Power))

# calculate mean across conditions across electrodes for each frequency/space condition
groupMean <- meanBySpace %>%
  group_by(Space, FreqBand) %>%
  summarise(GMean = mean(Power))

# normalize the observations
normBySpace <- inner_join(meanBySpace, elecMean)
normBySpace <- inner_join(normBySpace, groupMean) %>%
  mutate(NormPepisode = Power - EMean + GMean) %>%
  group_by(Space, FreqBand, TimePoint) %>%
  summarise(NormVar = var(NormPepisode)) %>%
  mutate(NormVarUnbias = NormVar * (length(timePointOrder) / (length(timePointOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(max(pepisode$Electrode)))) %>%
  select(Space, FreqBand, TimePoint, NormSEM)

# now summarise the original data and join with normalized SEM
meanBySpace <- meanBySpace %>%
  group_by(Space, FreqBand, TimePoint) %>%
  summarise(Pepisode = mean(Power)) %>%
  inner_join(normBySpace)

# plot the data
ggplot(meanBySpace, aes(x = TimePoint, y = Pepisode, ymin = Pepisode - NormSEM, ymax = Pepisode + NormSEM, color = Space)) +
  geom_point(size = 4) +
  geom_pointrange() +
  facet_wrap(~ FreqBand, scales = "free") +
  ggtitle("Pepisode by Space and Time Point")

ggsave('Figures/Pepisode_by_space_normSEM_11May2015.png')


####################################################
## summarize and plot the data by TIME condition ##
####################################################
meanByTime <- pepisode %>%
  mutate(Time = substr(as.character(Condition),3,4)) %>%
  select(-Condition) %>%
  group_by(Electrode, FreqBand,TimePoint)
timeOrder <- c('NT','FT')
meanByTime$Time <- factor(meanByTime$Time, levels = timeOrder)

# calculate mean across conditions for each electrode for each frequency
elecMean <- meanByTime %>%
  group_by(Time, Electrode, FreqBand) %>%
  summarise(EMean = mean(Power))

# calculate mean across conditions across electrodes for each frequency/time condition
groupMean <- meanByTime %>%
  group_by(Time, FreqBand) %>%
  summarise(GMean = mean(Power))

# normalize the observations
normByTime <- inner_join(meanByTime, elecMean)
normByTime <- inner_join(normByTime, groupMean) %>%
  mutate(NormPepisode = Power - EMean + GMean) %>%
  group_by(Time, FreqBand, TimePoint) %>%
  summarise(NormVar = var(NormPepisode)) %>%
  mutate(NormVarUnbias = NormVar * (length(timePointOrder) / (length(timePointOrder) - 1)),
         NormSEM = (sqrt(NormVarUnbias) / sqrt(max(pepisode$Electrode)))) %>%
  select(Time, FreqBand, TimePoint, NormSEM)

# now summarise the original data and join with normalized SEM
meanByTime <- meanByTime %>%
  group_by(Time, FreqBand, TimePoint) %>%
  summarise(Pepisode = mean(Power)) %>%
  inner_join(normByTime)

# plot the data
ggplot(meanByTime, aes(x = TimePoint, y = Pepisode, ymin = Pepisode - NormSEM, ymax = Pepisode + NormSEM, color = Time)) +
  geom_point(size = 4) +
  geom_pointrange() +
  facet_wrap(~ FreqBand, scales = "free") +
  ggtitle("Power by Time and Time Point")

ggsave('Figures/Pepisode_by_time_normSEM_11May2015.png')
