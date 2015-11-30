load('Rda/allCleanPowerData.Rda')

thisPower <- allPower %>%
  filter(ElectrodeID == 'UCDMC14_TeleporterB_LAD1')

summaryPower <- thisPower %>%
  ungroup() %>%
  group_by(Frequency, TimePoint) %>%
  summarise(MeanTele = mean(telePower),
            SEMTele  = sd(telePower) / sqrt(n()),
            MeanNav  = mean(navPower),
            SEMNav   = sd(navPower) / sqrt(n()))

pTele <- ggplot(summaryPower, aes(x = Frequency, y = MeanTele, ymin = MeanTele - SEMTele, ymax = MeanTele + SEMTele, color = TimePoint)) +
  geom_point() +
  geom_pointrange() +
  geom_line() +
  scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128)) +
  scale_y_log10() 

pNav <- ggplot(summaryPower, aes(x = Frequency, y = MeanNav, ymin = MeanNav - SEMNav, ymax = MeanNav + SEMNav, color = TimePoint)) +
  geom_point() +
  geom_pointrange() +
  geom_line() +
  scale_x_log10(breaks = c(2, 4, 8, 16, 32, 64, 128)) +
  scale_y_log10() 