function eegPath = getEegPath(subjectDir, subjectID, teleporter, electrode, suffix)
eegPath = [subjectDir subjectID '/Epoched Data/' subjectID '_' teleporter '_epoched_' electrode suffix];
end