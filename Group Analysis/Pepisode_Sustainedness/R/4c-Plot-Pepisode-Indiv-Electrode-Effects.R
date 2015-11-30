# Script Name:  4c-Plot-Pepisode-Indiv-Electrode-Effects.R
# Author:       Lindsay Vass
# Date:         14 September
# Purpose:      This script will load the data for pepisode during navigation and
#               teleportation and make scatterplots of the differences at each
#               time point for each electrode/frequency band.

library(dplyr)
library(ggplot2)
library(ggthemes)

load('Rda/allCleanData.Rda')
remove(navSustain, teleSustain, sessionInfo)


# Make scatterplot for pepisode -------------------------------------------

diffPepisode <- allPepisode %>%
  mutate(NavMinusTele = NavPepisode - TelePepisode) %>%
  group_by(ElectrodeID, TimePoint, FrequencyBand) %>%
  summarise(MeanNavMinusTele = mean(NavMinusTele),
            SEMNavMinusTele = sd(NavMinusTele) / sqrt(n()))

#colFun <- colorRampPalette(c("red", "orange", "black", "deepskyblue", "dodgerblue4"))
colFun <- colorRampPalette(c("dodgerblue4", "deepskyblue", "black", "orange", "red"))

diffPepisode$TimePoint <- plyr::mapvalues(diffPepisode$TimePoint, from = c("Pre1", "Tele", "Post1"), to = c("Pre", "Tele", "Post"))
diffPepisode$TimePoint <- factor(diffPepisode$TimePoint, levels = c("Pre", "Tele", "Post"))
diffPepisode$ElectrodeID <- factor(diffPepisode$ElectrodeID)

scatterPepisode <- diffPepisode %>%
  ggplot(aes(x = TimePoint, 
             y = MeanNavMinusTele, 
             ymin = MeanNavMinusTele - SEMNavMinusTele, 
             ymax = MeanNavMinusTele + SEMNavMinusTele, 
             group = ElectrodeID,
             colour = MeanNavMinusTele)) +
  geom_point(size = 5) +
  geom_pointrange() +
  geom_line() +
  scale_color_gradientn("Navigation - Teleporter", colours = colFun(5)) +
  theme_stata() +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(size = 30),
        legend.text = element_blank(),
        legend.title = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1.5),
        strip.background = element_rect(colour = "black", size = 0.75),
        panel.border = element_rect(colour = "black", size = 0.75, fill = NA),
        panel.grid.major.y = element_line(colour = "dimgray", linetype = "longdash")) +
  ylab(expression("Difference in Mean P"["Episode"])) +
  facet_grid(~ FrequencyBand) 
scatterPepisode
ggsave('Figures/SingleElectrodePepisode/Scatter_Nav_versus_Tele.pdf', useDingbats = FALSE, width = 16, height = 8)

# plot for each electrode
patientIDs <- data.frame(UCDMC = c('UCDMC13', 'UCDMC14', 'UCDMC15'),
                         P = c('P1', 'P2', 'P3'))
for (thisElec in 1:nlevels(diffPepisode$ElectrodeID)) {
  elec <- levels(diffPepisode$ElectrodeID)[thisElec]
  elecSplit <- strsplit(elec, '_')
  session <- ifelse(elecSplit[[1]][2] == 'TeleporterA', 'Session 1', 'Session 2')
  pTitle <- paste0(patientIDs$P[which(patientIDs$UCDMC == elecSplit[[1]][1])], ' ', session, ' ', elecSplit[[1]][3])
  
  p <- diffPepisode %>%
    filter(ElectrodeID == elec) %>%
    ggplot(aes(x = TimePoint, 
               y = MeanNavMinusTele, 
               ymin = MeanNavMinusTele - SEMNavMinusTele, 
               ymax = MeanNavMinusTele + SEMNavMinusTele, 
               group = ElectrodeID)) +
    geom_point(size = 5) +
    geom_pointrange() +
    geom_line() +
    theme_stata() +
    theme(plot.background = element_rect(fill = "white"),
          text = element_text(size = 30),
          legend.text = element_blank(),
          legend.title = element_text(size = 24),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1.5),
          strip.background = element_rect(colour = "black", size = 0.75),
          panel.border = element_rect(colour = "black", size = 0.75, fill = NA),
          panel.grid.major.y = element_line(colour = "dimgray", linetype = "longdash")) +
    ylab(expression("Difference in Mean P"["Episode"])) +
    facet_grid(~ FrequencyBand) +
    scale_y_continuous(limits = c(-0.2, 0.4)) +
    ggtitle(pTitle)
  p
  fileName <- paste0('Figures/SingleElectrodePepisode/EachElectrodeSeparate/', elec, '_Pepisode_Scatterplot.pdf')
  ggsave(fileName, useDingbats = FALSE, width = 16, height = 8)
}

numElec <- diffPepisode %>%
  ungroup() %>%
  select(ElectrodeID, FrequencyBand) %>%
  unique() %>%
  group_by(FrequencyBand) %>%
  summarise(Count = n())