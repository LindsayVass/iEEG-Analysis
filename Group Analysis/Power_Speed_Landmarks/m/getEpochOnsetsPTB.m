function [epochOnsets, removedOnsets] = getEpochOnsetsPTB(timeSyncPath, systemTimeStart, teleporterEpochsPath, samplingRate)
% epochOnsets = getEpochOnsetsPTB(timeSyncPath, systemTimeStart, teleporterEpochsPath)
%
% Purpose: Take the timing of epochs onsets in Unix ticks and convert them
%   to indices of the EEG data.
%
% INPUT:
%   timeSyncPath: path to the .mat file containing the parameters for the
%       linear regression relating ticks to EEG bins
%   systemTimeStart: vector of epoch onset times in Unix ticks; output from
%       "sampleBehavioralData.m"
%   teleporterEpochsPath: path to the mat file containing the onsets of the
%       teleportation epochs in unix ticks, to be excluded from analysis;
%       in the case of one EDF, there is one variable (epochsEDF)
%   samplingRate: sampling rate of the EEG data in Hz
%
% OUTPUT:
%   epochOnsets: vector of onset times in EEG indices
%   removedOnsets: indices of the onsets that were removed because they
%       overlap with teleportation
%
% Author: Lindsay Vass
% Date: 5 August 2015

load(timeSyncPath);

% load teleporter epochs info so we can exclude it later
load(teleporterEpochsPath);
teleEpochLength = table([1; 2], [1830; 2830], 'VariableNames', {'Index', 'LengthMs'});
epochTime       = table(eTime, epochsEDF, 'VariableNames', {'Index', 'OnsetBins'});
teleEpochInfo   = innerjoin(teleEpochLength, epochTime);
OffsetBins      = table(round(teleEpochInfo.OnsetBins + (teleEpochInfo.LengthMs * (1/1000) * samplingRate)), 'VariableNames', {'OffsetBins'});
teleEpochInfo   = [teleEpochInfo, OffsetBins];

epochOnsets = nan(size(systemTimeStart));

for thisInterval = 1:length(epochOnsets)
    epochOnsets(thisInterval) = round(systemTimeStart(thisInterval) * time_sync_regression(1) + time_sync_regression(2));
end

removedOnsets = [];
for thisTele = 1:length(teleEpochInfo.OnsetBins)
    thisOnset  = teleEpochInfo.OnsetBins(thisTele);
    thisOffset = teleEpochInfo.OffsetBins(thisTele);
    firstOnset = find(thisOnset > epochOnsets, 1, 'last' );
    lastOnset  = find(thisOffset < epochOnsets, 1 ) - 1;
    removedOnsets = cat(1, removedOnsets, [firstOnset:1:lastOnset]');
end

epochOnsets(removedOnsets) = [];

