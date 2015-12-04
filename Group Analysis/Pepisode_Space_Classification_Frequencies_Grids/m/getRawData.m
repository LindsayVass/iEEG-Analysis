function [allBinaryData, allEEGData, eegTimes] = getRawData(inputData, timesNTCalc, timesFTCalc, timesNT, timesFT, tol, subjectDir, frequencies)

% split electrodeID into subject, teleporter, electrode
subjectData = cell(length(inputData.ElectrodeID), 3);
electrodeID = inputData.ElectrodeID;
for thisRow = 1:length(electrodeID)
    subjectData(thisRow, :) = strsplit(electrodeID{thisRow}, '_');
end

allEEGData    = cell(length(inputData.ElectrodeID), 1);
allBinaryData = cell(length(inputData.ElectrodeID), 1);

for thisTrial = 1:length(inputData.TrialNumber)
    
    % load EEG data
    eegPath = [subjectDir subjectData{thisTrial, 1} '/Epoched Data/Pipeline/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_' subjectData{thisTrial, 3}(1:3) '_epoched.set'];
    EEG = pop_loadset(eegPath);
    
    % find indices corresponding to our time interval
    if strcmpi('NT', inputData.TrialTimeType(thisTrial)) == 1
        calcInds = find(EEG.times >= timesNTCalc(1) & EEG.times <= timesNTCalc(2));
        teleInds = find(EEG.times >= timesNT(1) & EEG.times <= timesNT(2));
    else
        calcInds = find(EEG.times >= timesFTCalc(1) & EEG.times <= timesFTCalc(2));
        teleInds = find(EEG.times >= timesFT(1) & EEG.times <= timesFT(2));
    end
    
    % Find the index of this channel in the EEG.data
    chanLabels = {EEG.chanlocs.labels};
    chanInd = find(strcmpi(subjectData{thisTrial, 3}, chanLabels));
    
    % select the eeg data for analysis
    eegData = EEG.data(chanInd, :, inputData.TrialNumber(thisTrial));
    allEEGData{thisTrial} = eegData(teleInds);
    
    % load the power distribution file
    powerDistSaveFile = [subjectDir subjectData{thisTrial, 1} '/Mat Files/Pepisode/Power Distributions/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_' subjectData{thisTrial, 3} '_power_distribution.mat'];
    load(powerDistSaveFile);
    
    % restrict to this frequency
    thisFrequency = inputData.Frequency(thisTrial);
    frequencyList = cell2mat({powerDistribution.frequency});
    freqInd       = find(frequencies == thisFrequency);
    
    [binaryMatrix, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistribution, frequencies, eegData, EEG.srate, 95, 3);
    
    %allBinaryData{thisTrial} = binaryMatrix(freqInd, teleInds);
    allBinaryData{thisTrial} = binaryMatrix(:, teleInds);
    
    if (abs(inputData.Pepisode(thisTrial)) - mean(binaryMatrix(freqInd, teleInds))) > tol
        error('Pepisode does not match expected value.')
    end
    
    eegTimes = EEG.times(teleInds);
    
end