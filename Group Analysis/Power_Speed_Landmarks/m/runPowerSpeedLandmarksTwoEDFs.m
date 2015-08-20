function [meanLogPowerVals, speedWhite, landmarkClean] = runPowerSpeedLandmarksTwoEDFs(behavioralPath, pulseTimingPaths, unepochedEEGPaths, frequencies, chanName, intervalMs, teleporterEpochsPath, eegSamplingRate)
% [meanLogPowerVals, speedWhite, landmarkClean] = runPowerSpeedLandmarksTwoEDFs(behavioralPath, pulseTimingPaths, unepochedEEGPaths, frequencies, chanName, intervalMs, teleporterEpochsPath, eegSamplingRate)
%
% Purpose: Perform an analysis to extract mean log(power) values, avatar
%   speed, and landmark richness for each segment of the navigation data.
%   This version of the analysis is to be run for sessions with two EDFs.
%   The txt file output by Unity is used to derive speed and landmark
%   richness. This function will return three regressors containing a value
%   for log(power), speed (whitened), and landmark richness, that
%   represents the mean for each segment of the data.
%
% INPUT:
%   behavioralPath: path to the txt file output from unity during active
%       navigation
%   pulseTimingPaths: two-value cell array containing the paths to each of
%       the pulse timing mat files (one for each EDF)
%   unepochedEEGPaths: two-value cell array containing the paths to each of
%       the EEGLAB datasets containing the unepoched data that has been
%       marked for artifacts (one for each EDF)
%   frequencies: vector of frequencies at which to extract power chanName:
%   string indicating the name of the electrode intervalMs: length of the
%   desired epoch in milliseconds (e.g., 200)
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
epochLabel = 3;

% Parse and sample text file
[systemTime, xPos, zPos] = parseBehavioralTxt(behavioralPath, 1);
[systemTimeStart, ~, landmarkTypeSamp, speedSamp] = sampleBehavioralData(systemTime, xPos, zPos, intervalMs);

% convert system time in ticks to EEG time in indices
[epochOnsets, EDF1inds, EDF2inds, ~, removedOnsets1, removedOnsets2] = getEpochOnsets(pulseTimingPaths, systemTimeStart, teleporterEpochsPath, eegSamplingRate);

% identify good epochs
epochOnsets1 = epochOnsets{1};
epochOnsets2 = epochOnsets{2};

unepochedEEG1Path = unepochedEEGPaths{1};
unepochedEEG2Path = unepochedEEGPaths{2};

[~, goodEpochs1] = makeMiniEpochs(epochOnsets1, epochLabel, epochInterval, unepochedEEG1Path);
[~, goodEpochs2] = makeMiniEpochs(epochOnsets2, epochLabel, epochInterval, unepochedEEG2Path);

% update regressors
[landmarkClean1, landmarkClean2] = cleanRegressor(landmarkTypeSamp, EDF1inds, EDF2inds, goodEpochs1, goodEpochs2, removedOnsets1, removedOnsets2);
[speedClean1, speedClean2]       = cleanRegressor(speedSamp, EDF1inds, EDF2inds, goodEpochs1, goodEpochs2, removedOnsets1, removedOnsets2);
speedWhite1 = prewhiten(speedClean1);
speedWhite2 = prewhiten(speedClean2);
epochOnsets1 = epochOnsets1(goodEpochs1);
epochOnsets2 = epochOnsets2(goodEpochs2);

% calculate power for one electrode
[logPowerVals1, EEG1] = calculateLogPower(unepochedEEG1Path, chanName, frequencies);
[logPowerVals2, EEG2] = calculateLogPower(unepochedEEG2Path, chanName, frequencies);

% get mean power for each epoch
meanLogPowerVals1 = extractMeanEpochLogPower(logPowerVals1, epochOnsets1, intervalMs, EEG1.srate);
meanLogPowerVals2 = extractMeanEpochLogPower(logPowerVals2, epochOnsets2, intervalMs, EEG2.srate);

% combine regressors
meanLogPowerVals = cat(2, meanLogPowerVals1, meanLogPowerVals2);
speedWhite       = cat(1, speedWhite1, speedWhite2);
landmarkClean    = cat(1, landmarkClean1, landmarkClean2);

end