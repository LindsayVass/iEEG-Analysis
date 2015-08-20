# Draw example raw trace with highlighted sections indicating valid oscillation

library(R.matlab)
library(ggplot2)
library(ggthemes)

oscData <- readMat('Figures/BAMM2015/BAMM_example_pepisode.mat')

rawTrace <- oscData$eegData
pepisodeBinary <- oscData$episodeVector
x <- seq_len(length(rawTrace))

rawTrace <- as.data.frame(t(rawTrace)) %>%
  mutate(x = x)
names(rawTrace) <- c('y', 'x')

episodes <- as.data.frame(t(pepisodeBinary) * 250) %>%
  mutate(x = x)
names(episodes) <- c('y', 'x')


blankP <- ggplot(rawTrace, aes(x = x, y = y)) + 
  geom_line() +
  theme_solid() +
  scale_y_continuous(limits = c(-250, 250))
ggsave('Figures/BAMM2015/exampleRawTrace.png')

episodeP <- ggplot(rawTrace, aes(x = x, y = y)) + 
  geom_polygon(data = episodes, aes(x = x, y = y), fill = "#d8161688") +
  geom_polygon(data = episodes, aes(x = x, y = -y + 0.15), fill = "#d8161688") +
  geom_line() +
  theme_solid() +
  scale_y_continuous(limits = c(-250, 250))
ggsave('Figures/BAMM2015/exampleRawTraceHighlighted.png')

pepisode <- sum(pepisodeBinary) / length(pepisodeBinary)
