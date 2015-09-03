

%% parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b/'))

% output from ../R/4c-Plot-Best-Electrode-Data.R
inputMat = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Space_Classification_Frequencies/mat/bestPepisodeTrials.mat';
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
subjectData = cell(length(bestData.ElectrodeID), 3);
electrodeID = bestData.ElectrodeID;
for thisRow = 1:length(electrodeID)
    subjectData(thisRow, :) = strsplit(electrodeID{thisRow}, '_');
end

allEEGData    = cell(length(bestData.ElectrodeID), 1);
allBinaryData = cell(length(bestData.ElectrodeID), 1);

for thisTrial = 1:length(bestData.TrialNumber)
    
    % load EEG data
    eegPath = [subjectDir subjectData{thisTrial, 1} '/Epoched Data/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_epoched_' subjectData{thisTrial, 3}(1:3) '_noSpikes_noWaves.set'];
    EEG = pop_loadset(eegPath);
    
    % find indices corresponding to our time interval
    if strcmpi('NT', bestData.TrialTimeType(thisTrial)) == 1
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
    eegData = EEG.data(chanInd, :, bestData.TrialNumber(thisTrial));
    allEEGData{thisTrial} = eegData(teleInds);
    
    % load the power distribution file
    powerDistSaveFile = [subjectDir subjectData{thisTrial, 1} '/Mat Files/Pepisode/Power Distributions/' subjectData{thisTrial, 1} '_' subjectData{thisTrial, 2} '_' subjectData{thisTrial, 3} '_power_distribution.mat'];
    load(powerDistSaveFile);
    
    % restrict to this frequency
    thisFrequency = bestData.Frequency(thisTrial);
    frequencyList = cell2mat({powerDistribution.frequency});
    %freqInd       = find(frequencyList == thisFrequency);
    freqInd       = find(frequencies == thisFrequency);
    
    [binaryMatrix, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistribution, frequencies, eegData, EEG.srate, 95, 3);
    
    %allBinaryData{thisTrial} = binaryMatrix(freqInd, teleInds);
    allBinaryData{thisTrial} = binaryMatrix(:, teleInds);
    
    if (abs(bestData.Pepisode(thisTrial)) - mean(binaryMatrix(freqInd, teleInds))) > tol
        error('Pepisode does not match expected value.')
    end
    
    eegTimes = EEG.times(teleInds);
    
end

save('../mat/RawDataBestPepisodeTrials_1-8Hz.mat', 'allEEGData', 'allBinaryData', 'eegTimes', 'frequencies')
