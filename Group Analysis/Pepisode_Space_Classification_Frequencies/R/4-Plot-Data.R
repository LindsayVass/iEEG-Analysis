# Script Name:  4-Plot-Data.R
# Author:       Lindsay Vass
# Date:         9 July 2015
# Purpose:      This script will plot the classification results from 
#               3-Analyze-Data.R

library(dplyr)
library(ggplot2)

theData <- meanClassification %>%
  group_by(Model) %>%
  summarise(SEM = sd(Accuracy) / sqrt(n()), Accuracy = mean(Accuracy))

p <- theData %>%
  ggplot(aes(x = Model, y = Accuracy, ymin = Accuracy - SEM, ymax = Accuracy + SEM)) +
  geom_bar(stat = "identity") +
  geom_errorbar() + 
  geom_hline(yintercept = 0.5, color = "red")