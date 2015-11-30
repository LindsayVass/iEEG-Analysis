# Load .mat files containing power values at each frequency and time point

library(R.matlab)
#dir.create('Rda')

filePath <- 'mat/UCDMC14_TeleporterB_LAD1.mat'
fileParts <- unlist(strsplit(substr(filePath, 5, nchar(filePath)), '_'))
thisData <- readMat(filePath)
powerVec <- as.vector(thisData$allPower)
freqs <- t(thisData$frequencies)
powerDf <- data.frame(Subject = fileParts[1], 
                      Session = fileParts[2],
                      Electrode = substr(fileParts[3], 1, 4),
                      Frequency = freqs,
                      TimePoint = rep(1:(length(powerVec) / length(freqs)), each = length(freqs)),
                      Power = powerVec)
saveFile <- paste0('Rda/', fileParts[1], '_', fileParts[2], '_', substr(fileParts[3], 1, 4), '.Rda')
save(powerDf, file = saveFile)
