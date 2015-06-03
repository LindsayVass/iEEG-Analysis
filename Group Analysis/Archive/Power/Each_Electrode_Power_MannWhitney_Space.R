# Script Name:  Each_Electrode_Power_MannWhitney_Space.R
# Author:       Lindsay Vass
# Date:         28 May 2015
# Purpose:      This script will test whether power significantly differs 
#               between spatial conditions (NS/FS) within each frequency band 
#               for each electrode. It will use Wilcoxon's signed rank test
#               to identify electrodes that show:
#                 1. NS > FS (Main Effect collapsed across timepoints)
#                 2. NS < FS (Main Effect collapsed across timepoint)
#               It will determine whether the number of electrodes showing a 
#               given effect is higher than expected by chance using a binomial
#               test. Finally, it will plot the frequencies of these effects as
#               a bar graph. NB Power is defined as z-scored log(Power), where 
#               the values are z-scored within each electrode and frequency band.

# Set up workspace and load data ------------------------------------------

library(dplyr)
library(reshape2)
library(permute)

setwd('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power/')
load("All_Subjects_Power_CLEAN_3sBuffer2015-05-27.Rda")


# Initialize output data frame and prepare input data frame ---------------

wilcoxonResults <- data.frame(ElectrodeID = NA,
                              FrequencyBand = NA,
                              Contrast = NA,
                              P = NA)
# Prepare input data
allPowerData <- allPowerData %>%
  mutate(TrialID = paste(EDF, TrialNumber, sep = "-")) %>% # Add a variable for TrialID since there can be a trial #1 for EDF1 and a trial #1 for EDF2
  filter(TimePoint == "Pre1" | TimePoint == "Tele" | TimePoint == "Post1") %>% # keep only timepoints of interest
  select(TrialSpaceType, TimePoint, Power:TrialID) 

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
      group_by(TrialID, TrialSpaceType) %>%
      summarise(Power = mean(Power))
    
    # Cast to wide format
    frequencyDataWide <- frequencyData %>%
      dcast(TrialID ~ TrialSpaceType, value.var = "Power")
    
    # Wilcoxon signed-rank tests
    NSGtFS  <- wilcox.test(frequencyDataWide$NS, frequencyDataWide$FS, alternative = "greater")
    NSLtFS  <- wilcox.test(frequencyDataWide$NS, frequencyDataWide$FS, alternative = "less")
    
    # Add to summary dataframe
    wilcoxonResults[thisRow, ]     <- c(levels(allPowerData$ElectrodeID)[thisElectrode], 
                                      levels(allPowerData$FrequencyBand)[thisFreqBand], 
                                      "NS > FS", 
                                      NSGtFS$p.value)
    wilcoxonResults[thisRow + 1, ] <- c(levels(allPowerData$ElectrodeID)[thisElectrode], 
                                      levels(allPowerData$FrequencyBand)[thisFreqBand], 
                                      "NS < FS", 
                                      NSLtFS$p.value)
    thisRow <- thisRow + 2
    
  } # end thisFreqBand
  
} # end thisElectrode


# Summarize significant wilcoxon results ----------------------------------

# Reorder frequency band levels
freqBandOrder <- c("Delta", "Theta", "Alpha", "Beta", "Gamma")
wilcoxonResults$FrequencyBand <- factor(wilcoxonResults$FrequencyBand, levels = freqBandOrder)

# Get the counts of the significant results
wilcoxonSigResults <- wilcoxonResults %>%
  filter(P < 0.05) %>%
  group_by(FrequencyBand, Contrast) %>%
  summarise(Count = n())

# Perform binomial test to determine whether counts are higher than expected by chance
binomP <- wilcoxonSigResults %>%
  rowwise() %>%
  do(BinomialTestP = binom.test(.$Count, 
                                nlevels(allPowerData$ElectrodeID), 
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
  geom_hline(aes(yintercept = 0.05*nlevels(allPowerData$ElectrodeID), 
                 color = "red")) +
  facet_grid(. ~ FrequencyBand) +
  geom_text(data = sigMarkers, label = "*", size = 18) +
  theme(strip.text.x = element_text(size = 18),
        plot.title = element_text(size = 24, vjust = 2),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18)) +
  ggtitle("Frequency of Significant Differences in Z-scored log(Power) by Frequency Band")
wilcoxonPlot

# Save the chart
today <- Sys.Date()
ggsave(filename = paste0("Figures/Each_Electrode_Power_MannWhitney_Space_Bar_", today, ".png"))


  
