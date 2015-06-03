# Script Name:  Each_Electrode_Epoched_Power_MannWhitney_Time.R
# Author:       Lindsay Vass
# Date:         29 May 2015
# Purpose:      This script will test whether power significantly differs 
#               between temporal conditions (NT/FT) within each frequency band 
#               for each electrode. It will use Wilcoxon's signed rank test
#               to identify electrodes that show:
#                 1. NT > FT (Main Effect collapsed across timepoints)
#                 2. NT < FT (Main Effect collapsed across timepoint)
#               It will determine whether the number of electrodes showing a 
#               given effect is higher than expected by chance using a binomial
#               test. Finally, it will plot the frequencies of these effects as
#               a bar graph. NB Power is defined as z-scored log(Power), where 
#               the values are z-scored within each electrode and frequency band.
#               This script differs from "Each_Electrode_Power_MannWhitney_Time"
#               in that the data is epoched prior to calculating power rather 
#               than calculating power before epoching.

# Set up workspace and load data ------------------------------------------

library(dplyr)
library(reshape2)
library(permute)

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
    
    csvFile <- paste0(expDir, subjectIDs[thisSubject], "/Mat Files/Power/Summary/", subjectIDs[thisSubject], "_", teleporters[thisTeleporter], "_epoched_power_summary_noSpikes_noWaves_3sBuffer.csv")
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
saveFile <- paste0("All_Subjects_Epoched_Power_CLEAN_3sBuffer_", today, ".Rda")
save(allPowerData, file = saveFile)


# Initialize output data frame and prepare input data frame ---------------

mannWhitneyResults <- data.frame(ElectrodeID = NA,
                              FrequencyBand = NA,
                              Contrast = NA,
                              P = NA)
thisRow <- 1

# Prepare input data
allPowerData <- allPowerData %>%
  filter(TimePoint == "Pre1" | TimePoint == "Tele" | TimePoint == "Post1") %>% # keep only timepoints of interest
  select(TrialNumber, TrialTimeType, TimePoint, Power:FrequencyBand) 

# Perform analysis for each electrode x frequency band --------------------

for (thisElectrode in 1:nlevels(allPowerData$ElectrodeID)){
  
  # Get all data for this electrode
  electrodeData <- allPowerData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElectrode]) %>%
    select(-ElectrodeID) 
  
  for (thisFreqBand in 1:nlevels(allPowerData$FrequencyBand)) {
    
    # Filter data for this frequency band
    frequencyData <- electrodeData %>%
      filter(FrequencyBand == levels(FrequencyBand)[thisFreqBand]) %>%
      select(-FrequencyBand) %>%
      group_by(TrialNumber, TrialTimeType) %>%
      summarise(Power = mean(Power))
    
    # Cast to wide format
    frequencyDataWide <- frequencyData %>%
      dcast(TrialNumber ~ TrialTimeType, value.var = "Power")
    
    # Wilcoxon signed-rank tests
    NTGtFT  <- wilcox.test(frequencyDataWide$NT, frequencyDataWide$FT, alternative = "greater")
    NTLtFT  <- wilcox.test(frequencyDataWide$NT, frequencyDataWide$FT, alternative = "less")
    
    # Add to summary dataframe
    mannWhitneyResults[thisRow, ]     <- c(levels(allPowerData$ElectrodeID)[thisElectrode], 
                                      levels(allPowerData$FrequencyBand)[thisFreqBand], 
                                      "NT > FT", 
                                      NTGtFT$p.value)
    mannWhitneyResults[thisRow + 1, ] <- c(levels(allPowerData$ElectrodeID)[thisElectrode], 
                                      levels(allPowerData$FrequencyBand)[thisFreqBand], 
                                      "NT < FT", 
                                      NTLtFT$p.value)
    thisRow <- thisRow + 2
    
  } # end thisFreqBand
  
} # end thisElectrode


# Summarize significant wilcoxon results ----------------------------------

# Reorder frequency band levels
freqBandOrder <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
mannWhitneyResults$FrequencyBand <- factor(mannWhitneyResults$FrequencyBand, levels = freqBandOrder)

# Get the counts of the significant results
mannWhitneySigResults <- mannWhitneyResults %>%
  filter(P < 0.05) %>%
  group_by(FrequencyBand, Contrast) %>%
  summarise(Count = n())

# Perform binomial test to determine whether counts are higher than expected by chance
binomP <- mannWhitneySigResults %>%
  rowwise() %>%
  do(BinomialTestP = binom.test(.$Count, 
                                nlevels(allPowerData$ElectrodeID), 
                                p = 0.05, 
                                alternative = "greater")$p.value)
mannWhitneySigResults <- cbind(mannWhitneySigResults, binomP) 

# Make a dataframe for the significance markers
sigMarkers <- mannWhitneySigResults %>%
  filter(BinomialTestP < 0.05) %>%
  mutate(Count = Count + 0.2) 

# Plot the bar charts
mannWhitneyPlot <- mannWhitneySigResults %>%
  ggplot(aes(x = Contrast, y = Count)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = 0.05*nlevels(allPowerData$ElectrodeID), 
                 color = "red")) +
  facet_grid(. ~ FrequencyBand) +
  geom_text(data = sigMarkers, label = "*", size = 18) +
  theme(strip.text.x = element_text(size = 18),
        plot.title = element_text(size = 24, vjust = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18)) +
  ggtitle("Frequency of Significant Differences in Z-scored log(Power) by Frequency Band")
mannWhitneyPlot

# Save the chart
today <- Sys.Date()
ggsave(filename = paste0("Figures/Each_Electrode_Epoched_Power_MannWhitney_Time_Bar_", today, ".png"))


# Follow-up Analyses ------------------------------------------------------

# The Mann-Whitney U tests indicated a significant difference in the Delta and
# Theta bands. To characterize this effect, we'll perform direct comparisons
# between NT & FT within each timepoint.

sigFreqBands <- factor(c("Delta","Theta"))

followUpResults <- data.frame(ElectrodeID = NA,
                              FrequencyBand = NA,
                              TimePoint = NA,
                              Contrast = NA,
                              P = NA)
thisRow <- 1

# Update levels of TimePoint since we've filtered many of them out
allPowerData$TimePoint <- factor(allPowerData$TimePoint)

for (thisElectrode in 1:nlevels(allPowerData$ElectrodeID)){
  
  # Get all data for this electrode
  electrodeData <- allPowerData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElectrode]) %>%
    select(-ElectrodeID) 
  
  
  for (thisFreqBand in 1:nlevels(sigFreqBands)) {
    
    # Filter data for this frequency band
    frequencyData <- electrodeData %>%
      filter(FrequencyBand == levels(sigFreqBands)[thisFreqBand]) %>%
      select(-FrequencyBand)
    
      
    for (thisTimePoint in 1:nlevels(allPowerData$TimePoint)) {
      
      # Filter data for this timepoint
      timePointData <- frequencyData %>%
        filter(TimePoint == levels(TimePoint)[thisTimePoint]) %>%
        select(-TimePoint) %>%
        group_by(TrialNumber, TrialTimeType) %>%
        summarise(Power = mean(Power))
      
      # Cast to wide format
      timePointDataWide <- timePointData %>%
        dcast(TrialNumber ~ TrialTimeType, value.var = "Power")
      
      # Mann Whitney U tests
      NTGtFT  <- wilcox.test(timePointDataWide$NT, timePointDataWide$FT, alternative = "greater")
      NTLtFT  <- wilcox.test(timePointDataWide$NT, timePointDataWide$FT, alternative = "less")
      
      # Add to summary dataframe
      followUpResults[thisRow, ]     <- c(levels(allPowerData$ElectrodeID)[thisElectrode],
                                          levels(allPowerData$FrequencyBand)[thisFreqBand],
                                          levels(allPowerData$TimePoint)[thisTimePoint],
                                          "NT > FT",
                                          NTGtFT$p.value)
      followUpResults[thisRow + 1, ] <- c(levels(allPowerData$ElectrodeID)[thisElectrode],
                                          levels(allPowerData$FrequencyBand)[thisFreqBand],
                                          levels(allPowerData$TimePoint)[thisTimePoint],
                                          "NT < FT",
                                          NTLtFT$p.value)
      thisRow <- thisRow + 2;
      
      
    } # end thisTimepoint

  } # end thisFreqBand
  
} # end thisElectrode
  


# Summarize significant follow-up analyses --------------------------------

# Reorder frequency band levels
followUpResults$FrequencyBand <- factor(followUpResults$FrequencyBand, levels = freqBandOrder)
followUpResults$TimePoint <- factor(followUpResults$TimePoint, levels = c("Pre1","Tele","Post1"))

# Get the counts of the significant results
followUpSigResults <- followUpResults %>%
  filter(P < 0.05) %>%
  group_by(FrequencyBand, TimePoint, Contrast) %>%
  summarise(Count = n())

# Perform binomial test to determine whether counts are higher than expected by chance
binomP <- followUpSigResults %>%
  rowwise() %>%
  do(BinomialTestP = binom.test(.$Count, 
                                nlevels(allPowerData$ElectrodeID), 
                                p = 0.05, 
                                alternative = "greater")$p.value)
followUpSigResults <- cbind(followUpSigResults, binomP) 

# Make a dataframe for the significance markers
sigMarkers <- followUpSigResults %>%
  filter(BinomialTestP < 0.05) %>%
  mutate(Count = Count + 0.2) 

# Plot the bar charts
followUpPlot <- followUpSigResults %>%
  ggplot(aes(x = Contrast, y = Count)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = 0.05*nlevels(allPowerData$ElectrodeID), 
                 color = "red")) +
  facet_grid(FrequencyBand ~ TimePoint) +
  geom_text(data = sigMarkers, label = "*", size = 18) +
  theme(strip.text = element_text(size = 18),
        plot.title = element_text(size = 24, vjust = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 15)) +
  ggtitle("Frequency of Significant Differences \nin Z-scored log(Power) by Frequency Band")
followUpPlot

# Save the chart
today <- Sys.Date()
ggsave(filename = paste0("Figures/Each_Electrode_Epoched_Power_MannWhitney_Time_FollowUp_By_TimePoint_Bar_", today, ".png"))

