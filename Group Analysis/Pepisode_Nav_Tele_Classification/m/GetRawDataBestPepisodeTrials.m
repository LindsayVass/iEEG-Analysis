

%% parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b/'))

% output from ../R/4-Select-Best-Trials.R
inputMat = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Nav_Tele_Classification/mat/bestPepisodeTrials.mat';
load(inputMat)


% segment of trial used to estimate pepisode (ms)
timesNTCalc = [-3000 4830];
timesFTCalc = [-3000 5830];

% segment of trial corresponding to teleportation (ms)
timesNT = [1 1830];
timesFT = [1 2830];

intervalsNT = timesNT(:,2) - timesNT(:, 1) + 1;
intervalsFT = timesFT(:,2) - timesFT(:, 1) + 1;
minInterval = min([intervalsNT; intervalsFT]);

durationThresh = 3; % duration for pepisode in cycles
minFrequency   = durationThresh / (minInterval / 1000);
maxFrequency   = 8;
frequencies    = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011
excludedFreqs  = frequencies(frequencies < minFrequency | frequencies > maxFrequency);
frequencies    = frequencies(frequencies >= minFrequency & frequencies <= maxFrequency);

tol = 0.005; % acceptable difference between previously measured pepisode and pepisode measured in this script


subjectDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';   

%% analysis

% split electrodeID into subject, teleporter, electrode
subjectData = cell(length(bestTele.ElectrodeID), 3);
electrodeID = bestTele.ElectrodeID;
for thisRow = 1:length(electrodeID)
    subjectData(thisRow, :) = strsplit(electrodeID{thisRow}, '_');
end

allTeleEEGData    = cell(length(bestTele.ElectrodeID), 1);
allTeleBinaryData = cell(length(bestTele.ElectrodeID), 1);

allNavEEGData    = cell(length(bestNav.ElectrodeID), 1);
allNavBinaryData = cell(length(bestNav.ElectrodeID), 1);

eegNtTimes = [];
eegFtTimes = [];

for thisTrial = 1:length(bestTele.TrialNumber)
    
    % load EEG data
    eegPath = [subjectDir subjectData{thisTrial, 1} '/Epoched Data/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_epoched_' subjectData{thisTrial, 3}(1:3) '_noSpikes_noWaves.set'];
    EEG = pop_loadset(eegPath);
    
    % find indices corresponding to our time interval
    if strcmpi('NT', bestTele.TrialTimeType(thisTrial)) == 1
        calcInds = find(EEG.times >= timesNTCalc(1) & EEG.times <= timesNTCalc(2));
        teleInds = find(EEG.times >= timesNT(1) & EEG.times <= timesNT(2));
        
        if isempty(eegNtTimes)
            eegNtTimes = EEG.times(teleInds);
        end
    else
        calcInds = find(EEG.times >= timesFTCalc(1) & EEG.times <= timesFTCalc(2));
        teleInds = find(EEG.times >= timesFT(1) & EEG.times <= timesFT(2));
        
        if isempty(eegFtTimes)
            eegFtTimes = EEG.times(teleInds);
        end
    end
    
    % Find the index of this channel in the EEG.data
    chanLabels = {EEG.chanlocs.labels};
    chanInd = find(strcmpi(subjectData{thisTrial, 3}, chanLabels));
    
    % select the eeg data for analysis
    eegData = EEG.data(chanInd, :, bestTele.TrialNumber(thisTrial));
    allTeleEEGData{thisTrial} = eegData(teleInds);
    
    % load the power distribution file
    powerDistSaveFile = [subjectDir subjectData{thisTrial, 1} '/Mat Files/Pepisode/Power Distributions/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_' subjectData{thisTrial, 3} '_power_distribution.mat'];
    load(powerDistSaveFile);
    
    % restrict to this frequency
    thisFrequency = bestTele.Frequency(thisTrial);
    frequencyList = cell2mat({powerDistribution.frequency});
    %freqInd       = find(frequencyList == thisFrequency);
    freqInd       = find(frequencies == thisFrequency);
    
    [binaryMatrix, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistribution, frequencies, eegData, EEG.srate, 95, 3);
    
    %allBinaryData{thisTrial} = binaryMatrix(freqInd, teleInds);
    allTeleBinaryData{thisTrial} = binaryMatrix(:, teleInds);
    
    if (abs(bestTele.Pepisode(thisTrial)) - mean(binaryMatrix(freqInd, teleInds))) > tol
        error('Pepisode does not match expected value.')
    end
    
end

for thisTrial = 1:length(bestNav.TrialNumber)
    
    % load EEG data
    eegPath = [subjectDir subjectData{thisTrial, 1} '/Epoched Data/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_epoched_' subjectData{thisTrial, 3}(1:3) '_navigation.set'];
    EEG = pop_loadset(eegPath);
    
    % find indices corresponding to our time interval
    if strcmpi('NT', bestNav.TrialTimeType(thisTrial)) == 1
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
    eegData = EEG.data(chanInd, :, bestNav.TrialNumber(thisTrial));
    allNavEEGData{thisTrial} = eegData(teleInds);
    
    % load the power distribution file
    powerDistSaveFile = [subjectDir subjectData{thisTrial, 1} '/Mat Files/Pepisode/Power Distributions/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_' subjectData{thisTrial, 3} '_power_distribution.mat'];
    load(powerDistSaveFile);
    
    % restrict to this frequency
    thisFrequency = bestNav.Frequency(thisTrial);
    frequencyList = cell2mat({powerDistribution.frequency});
    %freqInd       = find(frequencyList == thisFrequency);
    freqInd       = find(frequencies == thisFrequency);
    
    [binaryMatrix, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistribution, frequencies, eegData, EEG.srate, 95, 3);
    
    %allBinaryData{thisTrial} = binaryMatrix(freqInd, teleInds);
    allNavBinaryData{thisTrial} = binaryMatrix(:, teleInds);
    
    if (abs(bestNav.Pepisode(thisTrial)) - mean(binaryMatrix(freqInd, teleInds))) > tol
        error('Pepisode does not match expected value.')
    end
    
end

save('../mat/RawDataBestPepisodeTrials_1-8Hz.mat', 'allTeleEEGData', 'allTeleBinaryData',  'allNavEEGData', 'allNavBinaryData', 'eegNtTimes', 'eegFtTimes', 'frequencies')
