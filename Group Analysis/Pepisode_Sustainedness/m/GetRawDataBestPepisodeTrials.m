

%% parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b/'))

% output from ../R/4b-Select-Single-Trial-Data.R
%inputMat = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Sustainedness/mat/bestSingleTrials.mat';
inputMat = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Sustainedness/mat/testSingleTrials.mat';

load(inputMat)

% segments of trial used to estimate pepisode (ms)
timePointNames = {'Pre1' 'Tele' 'Post1'};

timesNT = [-1829 0; 1 1830; 1831 3660];
timesFT = [-2829 0; 1 2830; 2831 5660];

% frequencies to use
frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% Check that all frequencies are valid for pepisode
intervalsNT = timesNT(:,2) - timesNT(:, 1) + 1;
intervalsFT = timesFT(:,2) - timesFT(:, 1) + 1;
minInterval = min([intervalsNT; intervalsFT]);

durationThresh = 3; % duration for pepisode in cycles
minFrequency   = durationThresh / (minInterval / 1000);
maxFrequency   = 8;
excludedFreqs  = frequencies(frequencies < minFrequency | frequencies > maxFrequency);
frequencies    = frequencies(frequencies >= minFrequency & frequencies <= maxFrequency);

tol = 0.005; % acceptable difference between previously measured pepisode and pepisode measured in this script


subjectDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';

%% analysis

% split electrodeID into subject, teleporter, electrode
subjectData = cell(length(manualBestTrials.ElectrodeID), 3);
electrodeID = manualBestTrials.ElectrodeID;
for thisRow = 1:length(electrodeID)
    subjectData(thisRow, :) = strsplit(electrodeID{thisRow}, '_');
end



conditions = {'nav', 'tele'};
allEEGData    = cell(length(manualBestTrials.ElectrodeID), length(conditions));
allBinaryData = cell(length(manualBestTrials.ElectrodeID), length(conditions));
allEEGTimes   = cell(length(manualBestTrials.ElectrodeID), length(conditions));

for thisObs = 1:length(manualBestTrials.RealTrialNumber)
    for navTele = 1:length(conditions)
        
        % load EEG data
        if strcmpi(conditions{navTele}, 'nav') == 1
            suffix = '_navigation.set';
        else
            suffix = '_noSpikes_noWaves.set';
        end
        
        eegPath = getEegPath(subjectDir, subjectData{thisObs, 1}, subjectData{thisObs, 2}, subjectData{thisObs, 3}(1:3), suffix);
        EEG = pop_loadset(eegPath);
        
        % find indices corresponding to our time interval
        timesInd = find(strcmpi(manualBestTrials.TimePoint{thisObs}, timePointNames) == 1);
        if strcmpi('NT', manualBestTrials.TrialTimeType(thisObs)) == 1
            eegInds = find(EEG.times >= timesNT(timesInd, 1) & EEG.times <= timesNT(timesInd, 2));
        else
            eegInds = find(EEG.times >= timesFT(timesInd, 1) & EEG.times <= timesFT(timesInd, 2));
        end
        
        % Find the index of this channel in the EEG.data
        chanLabels = {EEG.chanlocs.labels};
        chanInd = find(strcmpi(subjectData{thisObs, 3}, chanLabels));
        
        % select the eeg data for analysis
        if strcmpi(conditions{navTele}, 'nav') == 1
            eegData = EEG.data(chanInd, eegInds, manualBestTrials.NavTrialNumber(thisObs));
        else
            eegData = EEG.data(chanInd, eegInds, manualBestTrials.TeleTrialNumber(thisObs));
        end
        allEEGData{thisObs, navTele} = eegData;
        
        % load the power distribution file
        powerDistSaveFile = [subjectDir subjectData{thisObs, 1} '/Mat Files/Pepisode/Power Distributions/' subjectData{thisObs, 1} '_' subjectData{thisObs, 2} '_' subjectData{thisObs, 3} '_power_distribution.mat'];
        load(powerDistSaveFile);
        
        % restrict to this frequency
        thisFrequency = manualBestTrials.Frequency(thisObs);
        freqInd       = find(abs(frequencies - thisFrequency) < tol);
        
        [binaryMatrix, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistribution, frequencies, eegData, EEG.srate, 95, 3);
        
        allBinaryData{thisObs, navTele} = binaryMatrix(freqInd, :);
        
        if strcmpi(conditions{navTele}, 'nav') == 1
            targetPepisode = manualBestTrials.NavPepisode(thisObs);
        else
            targetPepisode = manualBestTrials.TelePepisode(thisObs);
        end
        
        if (abs(targetPepisode - mean(binaryMatrix(freqInd, :)))) > tol
            error('Pepisode does not match expected value.')
        end
        
        eegTimes = EEG.times(eegInds);
    end
end

%save('../mat/bestSingleTrialsRawData.mat', 'allEEGData', 'allBinaryData', 'allEEGTimes')
