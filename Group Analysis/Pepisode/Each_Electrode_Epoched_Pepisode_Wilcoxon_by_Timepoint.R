# Script Name:  Each_Electrode_Epoched_Pepisode_Wilcoxon_by_Timepoint.R
# Author:       Lindsay Vass
# Date:         2 June 2015
# Purpose:      This script will test whether pepisode significantly differs 
#               between timepoints (Pre1, Tele, Post1) within each frequency
#               band for each electrode. It will use Wilcoxon signed rank tests
#               to identify electrodes that show:
#                 1. Pre > Tele
#                 2. Pre < Tele
#                 3. Tele > Post
#                 4. Tele < Post
#               It will determine whether the number of electrodes showing a 
#               given effect is higher than expected by chance using a binomial
#               test. Finally, it will plot the frequencies of these effects as
#               a bar graph. Data were epoched prior to calculating pepisode.

# Set up workspace ------------------------------------------------------------

library(dplyr)
library(reshape2)
library(permute)
library(ggplot2)
library(tidyr)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/')


# Load in the data for each session for each subject ----------------------
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
    
    csvFile <- paste0(expDir, subjectIDs[thisSubject], "/Mat Files/Pepisode/Summary/", subjectIDs[thisSubject], "_", teleporters[thisTeleporter], "_epoched_pepisode_summary_noSpikes_noWaves_3sBuffer.csv")
    tempData <- read.csv(csvFile, header = TRUE)
    
    allPepisodeData <- rbind(allPepisodeData, tempData)
    
  } # thisTeleporter
  
} # thisSubject

## Put our conditions in order
timePointOrder <- c('Pre3', 'Pre2', 'Pre1', 'Tele', 'Post1', 'Post2', 'Post3')
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
freqBandBreaks <- c(0, 8, 12, 30, 182)
freqBandNames  <- c("Delta-Theta","Alpha","Beta","Gamma")

allPepisodeData <- allPepisodeData %>%
  mutate(FrequencyBand = cut(Frequency, freqBandBreaks, labels = freqBandNames))

# Write the data out for stats later
today <- Sys.Date()
saveFile <- paste0("All_Subjects_Epoched_Pepisode_CLEAN_3sBuffer_", today, ".Rda")
save(allPepisodeData, file = saveFile)


# Initialize output data frame and prepare input data frame ---------------

wilcoxonResults <- data.frame(ElectrodeID = NA,
                              FrequencyBand = NA,
                              Contrast = NA,
                              P = NA)
thisRow <- 1

# Prepare input data
allPepisodeData <- allPepisodeData %>%
  filter(TimePoint == "Pre1" | TimePoint == "Tele" | TimePoint == "Post1") %>% # keep only timepoints of interest
  select(TrialNumber, TimePoint, Pepisode:FrequencyBand) 

# Perform analysis for each electrode x frequency band --------------------

for (thisElectrode in 1:nlevels(allPepisodeData$ElectrodeID)){
  
  # Get all data for this electrode
  electrodeData <- allPepisodeData %>%
    filter(ElectrodeID == levels(ElectrodeID)[thisElectrode]) %>%
    select(-ElectrodeID) 
  
  for (thisFreqBand in 1:nlevels(allPepisodeData$FrequencyBand)) {
    
    # Filter data for this frequency band
    frequencyData <- electrodeData %>%
      filter(FrequencyBand == levels(FrequencyBand)[thisFreqBand]) %>%
      select(-FrequencyBand) %>%
      group_by(TrialNumber, TimePoint) %>%
      summarise(Pepisode = mean(Pepisode))
    
    # Cast to wide format
    frequencyDataWide <- frequencyData %>%
      dcast(TrialNumber ~ TimePoint, value.var = "Pepisode")
    
    # Wilcoxon signed-rank tests
    preGtTele  <- wilcox.test(frequencyDataWide$Pre1, frequencyDataWide$Tele, alternative = "greater", paired=TRUE)
    preLtTele  <- wilcox.test(frequencyDataWide$Pre1, frequencyDataWide$Tele, alternative = "less", paired=TRUE)
    teleGtPost <- wilcox.test(frequencyDataWide$Tele, frequencyDataWide$Post1, alternative = "greater", paired = TRUE)
    teleLtPost <- wilcox.test(frequencyDataWide$Tele, frequencyDataWide$Post1, alternative = "less", paired = TRUE)
    
    # Add to summary dataframe
    wilcoxonResults[thisRow, ]     <- c(levels(allPepisodeData$ElectrodeID)[thisElectrode], 
                                      levels(allPepisodeData$FrequencyBand)[thisFreqBand], 
                                      "Pre > Tele", 
                                      preGtTele$p.value)
    wilcoxonResults[thisRow + 1, ] <- c(levels(allPepisodeData$ElectrodeID)[thisElectrode], 
                                      levels(allPepisodeData$FrequencyBand)[thisFreqBand], 
                                      "Pre < Tele", 
                                      preLtTele$p.value)
    wilcoxonResults[thisRow + 2, ] <- c(levels(allPepisodeData$ElectrodeID)[thisElectrode], 
                                      levels(allPepisodeData$FrequencyBand)[thisFreqBand], 
                                      "Tele > Post", 
                                      teleGtPost$p.value)
    wilcoxonResults[thisRow + 3, ] <- c(levels(allPepisodeData$ElectrodeID)[thisElectrode], 
                                      levels(allPepisodeData$FrequencyBand)[thisFreqBand], 
                                      "Tele < Post", 
                                      teleLtPost$p.value)
    thisRow <- thisRow + 4
    
  } # end thisFreqBand
  
} # end thisElectrode


# Summarize significant wilcoxon results ----------------------------------

# Reorder frequency band levels
freqBandOrder <- c("Delta-Theta", "Alpha", "Beta", "Gamma")
wilcoxonResults$FrequencyBand <- factor(wilcoxonResults$FrequencyBand, levels = freqBandOrder)

# Get the counts of the significant results
wilcoxonSigResults <- wilcoxonResults %>%
  group_by(FrequencyBand, Contrast) %>%
  filter(P < 0.05) %>%
  summarise(Count = n())

# Fill in missing frequency band/contrasts with zeros
allContrasts <- wilcoxonResults %>%
  group_by(FrequencyBand, Contrast) %>%
  summarise(Count = 0) %>%
  anti_join(wilcoxonSigResults, by = c("FrequencyBand", "Contrast"))
wilcoxonSigResults <- rbind(wilcoxonSigResults, allContrasts)


# Perform binomial test to determine whether counts are higher than expected by chance
binomP <- wilcoxonSigResults %>%
  rowwise() %>%
  do(BinomialTestP = binom.test(.$Count, 
                                nlevels(allPepisodeData$ElectrodeID), 
                                p = 0.05, 
                                alternative = "greater")$p.value)
wilcoxonSigResults <- cbind(wilcoxonSigResults, binomP) 

# Make a dataframe for the significance markers
sigMarkers <- wilcoxonSigResults %>%
  filter(BinomialTestP < 0.05) %>%
  mutate(Count = Count + 0.2) 

# Plot the bar charts
wilcoxonPlot <- wilcoxonSigResults %>%
  ggplot(aes(x = Contrast, y = Count)) +
  geom_bar(stat = "identity") +
  geom_hline(aes(yintercept = 0.05*nlevels(allPepisodeData$ElectrodeID), 
                 color = "red")) +
  facet_grid(. ~ FrequencyBand) +
  geom_text(data = sigMarkers, label = "*", size = 18) +
  theme(strip.text.x = element_text(size = 18),
        plot.title = element_text(size = 24, vjust = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18)) +
  ggtitle("Frequency of Significant Differences in Pepisode by Frequency Band")
wilcoxonPlot

# Save the chart
today <- Sys.Date()
ggsave(filename = paste0("Figures/Each_Electrode_Epoched_Pepisode_Wilcoxon_by_Timepoint_Bar_", today, ".png"))

# Plot data for each electrode separately ---------------------------------
for (thisFreqBand in 1:nlevels(allPepisodeData$FrequencyBand)) {
  electrodeWisePlot <- allPepisodeData %>%
    filter(FrequencyBand == levels(allPepisodeData$FrequencyBand)[thisFreqBand]) %>%
    group_by(TimePoint, ElectrodeID) %>%
    mutate(SEM = sd(Pepisode) / sqrt(n() - 1)) %>%
    summarise(Pepisode = mean(Pepisode), SEM = mean(SEM)) %>%
    ggplot(aes(x = TimePoint, y = Pepisode, ymin = Pepisode - SEM, ymax = Pepisode + SEM)) +
    geom_point(size = 4) +
    geom_pointrange() +
    facet_wrap(~ ElectrodeID, scales = "free") +
    ggtitle(paste0("Z-Scored ", levels(allPepisodeData$FrequencyBand)[thisFreqBand], " Pepisode"))
  
  electrodeWisePlot
  
  # Save the chart
  today <- Sys.Date()
  ggsave(filename = paste0("Figures/Each_Electrode_Epoched_", levels(allPepisodeData$FrequencyBand)[thisFreqBand], "_Pepisode_by_Timepoint_", today, ".png"))
  
}

