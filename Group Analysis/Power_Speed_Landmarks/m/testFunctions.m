% %% Test 1 EDF (UCDMC13 & UCDMC14)
% clear all; close all; clc;
% 
% % parameters
% intervalMs = 200;
% epochLabel = 3;
% epochInterval = [0 0.2];
% 
% % Parse and sample text file
% behavioralPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Behavioral Data/TeleporterA/s2_patientTeleporterData.txt';
% [systemTime, xPos, zPos] = parseBehavioralTxt(behavioralPath);
% [systemTimeStart, systemTimeEnd, landmarkTypeSamp, speedSamp] = sampleBehavioralData(systemTime, xPos, zPos, intervalMs);
% 
% assert((length(systemTimeStart) == length(systemTimeEnd)), 'Sampling returned regressors of different lengths. (systemTimeStart & systemTimeEnd)')
% assert((length(systemTimeEnd) == length(landmarkTypeSamp)), 'Sampling returned regressors of different lengths. (systemTimeEnd & landmarkTypeSamp)')
% assert((length(landmarkTypeSamp) == length(speedSamp)), 'Sampling returned regressors of different lengths. (landmarkTypeSamp & speedSamp)')
% assert(length(systemTimeStart) < length(systemTime), 'Sampling did not reduce the number of values in systemTimeStart.')
% 
% % convert system time in ticks to EEG time in indices
% timeSyncPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Mat Files/UCDMC14_TeleporterA_time_sync.mat';
% teleporterEpochsPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Mat Files/UCDMC14_TeleporterA_Epochs_Entry.mat';
% eegSamplingRate = 512;
% [epochOnsets, removedOnsets] = getEpochOnsetsPTB(timeSyncPath, systemTimeStart, teleporterEpochsPath, eegSamplingRate);
% 
% assert(length(epochOnsets) + length(removedOnsets) == length(systemTimeStart), 'Converting ticks to EEG indices changed the number of time points.')
% systemTimeStart(removedOnsets) = [];
% diffTicks = diff(systemTimeStart);
% diffEEG   = diff(epochOnsets);
% diffCorr  = corrcoef(diffTicks, diffEEG) > 0.9;
% assert(sum(diffCorr(:)) == 4, 'Tick intervals are not sufficiently correlated with EEG intervals.')
% 
% % identify good epochs
% unepochedEEGPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/PreProcessing Intermediates/UCDMC14_TeleporterA_unepoched_LAD_marked.set';
% [epochedEEG, goodEpochs] = makeMiniEpochs(epochOnsets, epochLabel, epochInterval, unepochedEEGPath);
% 
% assert(length(goodEpochs) <= length(epochOnsets), 'More epochs in clean data than original data.')
% 
% % update regressors
% landmarkClean = landmarkTypeSamp(goodEpochs);
% speedClean    = speedSamp(goodEpochs);
% speedWhite    = prewhiten(speedClean);
% epochOnsets   = epochOnsets(goodEpochs);
% 
% % calculate power for one electrode
% frequencies  = logspace(log(1)/log(10),log(181)/log(10),31); 
% chanName     = 'LAD1';
% logPowerVals = calculateLogPower(unepochedEEGPath, chanName, frequencies);
% EEG = pop_loadset(unepochedEEGPath);
% 
% assert(size(logPowerVals, 1) == length(frequencies), 'Power calculation did not return the correct number of frequencies.')
% assert(size(logPowerVals, 2) == size(EEG.data, 2), 'Power calculation did not return the correct number of time points.')
% 
% % get mean power for each epoch
% meanLogPowerVals = extractMeanEpochLogPower(logPowerVals, epochOnsets, intervalMs, EEG.srate);
% 
% assert(size(meanLogPowerVals, 1) == length(frequencies), 'Mean epoch power did not return the correct number of frequencies.')
% assert(size(meanLogPowerVals, 2) == length(epochOnsets), 'Mean epoch power did not return the correct number of epochs.')
% epochLengthBins = round(intervalMs * (1/1000) * EEG.srate) - 1; % subtract 1 because epochOnset will be the first bin
% oneFreqOneEpoch = mean(logPowerVals(5, epochOnsets(5):epochOnsets(5) + epochLengthBins));
% assert(oneFreqOneEpoch == meanLogPowerVals(5,5), 'Mean epoch power does not match expected value.')


%% Test 2 EDFs (UCDMC15)
clear all; close all; 

% parameters
intervalMs = 200;
epochLabel = 3;
epochInterval = [0 0.2];

% Parse and sample text file
behavioralPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Behavioral Data/TeleporterA/s3_FindStore_TeleporterA_FIXED.txt';
[systemTime, xPos, zPos] = parseBehavioralTxt(behavioralPath, 1);
[systemTimeStart, systemTimeEnd, landmarkTypeSamp, speedSamp] = sampleBehavioralData(systemTime, xPos, zPos, intervalMs);
assert((length(systemTimeStart) == length(systemTimeEnd)), 'Sampling returned regressors of different lengths. (systemTimeStart & systemTimeEnd)')
assert((length(systemTimeEnd) == length(landmarkTypeSamp)), 'Sampling returned regressors of different lengths. (systemTimeEnd & landmarkTypeSamp)')
assert((length(landmarkTypeSamp) == length(speedSamp)), 'Sampling returned regressors of different lengths. (landmarkTypeSamp & speedSamp)')
assert(length(systemTimeStart) < length(systemTime), 'Sampling did not reduce the number of values in systemTimeStart.')


% convert system time in ticks to EEG time in indices
pulseTimingPaths = {'/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Mat Files/UCDMC15_TeleporterA_EDF1_pulse_timing.mat', ...
    '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Mat Files/UCDMC15_TeleporterA_EDF2_pulse_timing.mat'};
teleporterEpochsPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Mat Files/UCDMC15_TeleporterA_Epochs_Entry.mat';
eegSamplingRate = 512;
[epochOnsets, EDF1inds, EDF2inds, missingInds, removedOnsets1, removedOnsets2] = getEpochOnsets(pulseTimingPaths, systemTimeStart, teleporterEpochsPath, eegSamplingRate);

assert(size(epochOnsets{1}, 2) + size(epochOnsets{2}, 2) + length(missingInds) + length(removedOnsets1) + length(removedOnsets2) == length(systemTimeStart), 'Converting ticks to EEG indices changed the number of time points.')
newTicks1  = systemTimeStart(EDF1inds);
newTicks1(removedOnsets1) = [];
newTicks2  = systemTimeStart(EDF2inds);
newTicks2(removedOnsets2) = [];
diffTicks1 = diff(newTicks1);
diffTicks2 = diff(newTicks2);
diffTicks  = cat(1, diffTicks1, diffTicks2);
diffEEG1  = diff(epochOnsets{1})';
diffEEG2  = diff(epochOnsets{2})';
diffEEG   = cat(1, diffEEG1, diffEEG2);
diffCorr  = corrcoef(diffTicks, diffEEG) > 0.9; % this one is not as highly correlated as the PTB method, probably because that's using the same linear regression, whereas this one calculates a new fit for each data point based on local pulses
assert(sum(diffCorr(:)) == 4, 'Tick intervals are not sufficiently correlated with EEG intervals.')

% identify good epochs
unepochedEEG1Path = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/PreProcessing Intermediates/UCDMC15_TeleporterA_EDF1_unepoched_LAD_marked.set';
unepochedEEG2Path = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/PreProcessing Intermediates/UCDMC15_TeleporterA_EDF2_unepoched_LAD_marked.set';

epochOnsets1 = epochOnsets{1};
epochOnsets2 = epochOnsets{2};

[epochedEEG1, goodEpochs1] = makeMiniEpochs(epochOnsets1, epochLabel, epochInterval, unepochedEEG1Path);
[epochedEEG2, goodEpochs2] = makeMiniEpochs(epochOnsets2, epochLabel, epochInterval, unepochedEEG2Path);

assert(length(goodEpochs1) <= length(epochOnsets1), 'More epochs in EDF1 clean data than original data.')
assert(length(goodEpochs2) <= length(epochOnsets2), 'More epochs in EDF2 clean data than original data.')

% update regressors
[landmarkClean1, landmarkClean2] = cleanRegressor(landmarkTypeSamp, EDF1inds, EDF2inds, goodEpochs1, goodEpochs2, removedOnsets1, removedOnsets2);
[speedClean1, speedClean2]       = cleanRegressor(speedSamp, EDF1inds, EDF2inds, goodEpochs1, goodEpochs2, removedOnsets1, removedOnsets2);
speedWhite1 = prewhiten(speedClean1);
speedWhite2 = prewhiten(speedClean2);
epochOnsets1 = epochOnsets1(goodEpochs1);
epochOnsets2 = epochOnsets2(goodEpochs2);

assert(length(epochOnsets1) == length(goodEpochs1), 'Number of good epochs for EDF1 is incorrect.')
assert(length(epochOnsets2) == length(goodEpochs2), 'Number of good epochs for EDF2 is incorrect.')

% calculate power for one electrode
frequencies  = logspace(log(1)/log(10),log(181)/log(10),31); 
chanName     = 'LAD1';
logPowerVals1 = calculateLogPower(unepochedEEG1Path, chanName, frequencies);
logPowerVals2 = calculateLogPower(unepochedEEG2Path, chanName, frequencies);
EEG1 = pop_loadset(unepochedEEG1Path);
EEG2 = pop_loadset(unepochedEEG2Path);

assert(size(logPowerVals1, 1) == length(frequencies), 'Power calculation did not return the correct number of frequencies.')
assert(size(logPowerVals1, 2) + size(logPowerVals2, 2) == size(EEG1.data, 2) + size(EEG2.data, 2), 'Power calculation did not return the correct number of time points.')

% get mean power for each epoch
meanLogPowerVals = extractMeanEpochLogPower(logPowerVals1, epochOnsets1, intervalMs, EEG1.srate);

assert(size(meanLogPowerVals, 1) == length(frequencies), 'Mean epoch power did not return the correct number of frequencies.')
assert(size(meanLogPowerVals, 2) == length(epochOnsets1), 'Mean epoch power did not return the correct number of epochs.')
epochLengthBins = round(intervalMs * (1/1000) * EEG1.srate) - 1; % subtract 1 because epochOnset will be the first bin
oneFreqOneEpoch = mean(logPowerVals1(5, epochOnsets1(5):epochOnsets1(5) + epochLengthBins));
assert(oneFreqOneEpoch == meanLogPowerVals(5,5), 'Mean epoch power does not match expected value.')

