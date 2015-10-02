# Script Name:  3-Plot-Data.R
# Author:       Lindsay Vass
# Date:         1 October 2015
# Purpose:      This script will plot the mean speed for each electrode for Pre
#               vs Post.

library(ggplot2)

load('Rda/allAnalyzedData.Rda')
#dir.create('Figures')

summarySpeedData$Interval <- factor(summarySpeedData$Interval, levels = c("Pre", "Post"))

speedBar <- summarySpeedData %>%
  ggplot(aes(x = Interval, y = MeanSpeed, ymin = MeanSpeed - SEMSpeed, ymax = MeanSpeed + SEMSpeed)) +
  geom_bar(stat = "identity") +
  geom_errorbar() +
  facet_wrap(~ElectrodeID) +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank())
speedBar

ggsave('Figures/Bar_Speed_Pre_Post.png')