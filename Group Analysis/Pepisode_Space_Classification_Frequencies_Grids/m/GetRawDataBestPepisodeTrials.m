

%% parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b/'))

% output from ../R/5-Select-Raw-Data.R
inputMat = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Space_Classification_Frequencies_Grids/mat/pepisodeTrials.mat';
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

[allBinaryData, allEEGData, eegTimes] = getRawData(occipData, timesNTCalc, timesFTCalc, timesNT, timesFT, tol, subjectDir, frequencies);
save('../mat/OccipRawDataPepisodeTrials_1-8Hz.mat', 'allEEGData', 'allBinaryData', 'eegTimes', 'frequencies')

[allBinaryData, allEEGData, eegTimes] = getRawData(parieData, timesNTCalc, timesFTCalc, timesNT, timesFT, tol, subjectDir, frequencies);
save('../mat/ParieRawDataPepisodeTrials_1-8Hz.mat', 'allEEGData', 'allBinaryData', 'eegTimes', 'frequencies')

