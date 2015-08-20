function [meanLogPowerVals, speedWhite, landmarkClean] = runPowerSpeedLandmarksOneEDF(behavioralPath, timeSyncPath, unepochedEEGPath, frequencies, chanName, intervalMs, teleporterEpochsPath, eegSamplingRate)
% [meanLogPowerVals, speedWhite, landmarkClean] = runPowerSpeedLandmarksOneEDF(behavioralPath, timeSyncPath, unepochedEEGPath, frequencies, chanName, intervalMs, teleporterEpochsPath, eegSamplingRate)
%
% Purpose: Perform an analysis to extract mean log(power) values, avatar
%   speed, and landmark richness for each segment of the navigation data.
%   This version of the analysis is to be run for sessions with one EDF.
%   The txt file output by Unity is used to derive speed and landmark
%   richness. This function will return three regressors containing a value
%   for log(power), speed (whitened), and landmark richness, that
%   represents the mean for each segment of the data.
%
% INPUT:
%   behavioralPath: path to the txt file output from unity during active
%       navigation
%   timeSyncPath: path to the mat file that contains the result of the
%       linear regression to map unix ticks to EEG indices
%   unepochedEEGPath: path to the EEGLAB dataset containing the unepoched
%       data that has been marked for artifacts
%   frequencies: vector of frequencies at which to extract power
%   chanName: string indicating the name of the electrode
%   intervalMs: length of the desired epoch in milliseconds (e.g., 200)
%   teleporterEpochsPath: path to the mat file containing the onsets of the
%       teleportation epochs in unix ticks, to be excluded from analysis;
%       in the case of one EDF, there is one variable (epochsEDF); in the
%       case of two EDFs, there are two variables (epochsEDF1 & epochsEDF2)
%   eegSamplingRate: sampling rate of the EEG data in Hz
%
% OUTPUT:
%   meanLogPowerVals: frequencies x epochs matrix containing the mean
%       log(Power) for each epoch of the EEG
%   speedWhite: prewhitened vector containing the mean speed for each epoch
%       of the EEG
%   landmarkClean: cell array containing the modal landmark richness for
%       each epoch of the EEG
%
% Author: Lindsay Vass
% Date: 6 August 2015

epochInterval = [0 intervalMs / 1000];

% Parse and sample text file
[systemTime, xPos, zPos] = parseBehavioralTxt(behavioralPath);
[systemTimeStart, ~, landmarkTypeSamp, speedSamp] = sampleBehavioralData(systemTime, xPos, zPos, intervalMs);

% convert system time in ticks to EEG time in indices
[epochOnsets, removedOnsets] = getEpochOnsetsPTB(timeSyncPath, systemTimeStart, teleporterEpochsPath, eegSamplingRate);

% remove epochs from our regressors that overlapped with teleportation
landmarkTypeSamp(removedOnsets) = [];
speedSamp(removedOnsets)        = [];

% identify good epochs
epochLabel = 3;
[~, goodEpochs] = makeMiniEpochs(epochOnsets, epochLabel, epochInterval, unepochedEEGPath);

% update regressors
landmarkClean = landmarkTypeSamp(goodEpochs);
speedClean    = speedSamp(goodEpochs);
speedWhite    = prewhiten(speedClean);
epochOnsets   = epochOnsets(goodEpochs);

% calculate power for one electrode
[logPowerVals, EEG] = calculateLogPower(unepochedEEGPath, chanName, frequencies);

% get mean power for each epoch
meanLogPowerVals = extractMeanEpochLogPower(logPowerVals, epochOnsets, intervalMs, EEG.srate);

end