function [epochOnsets, EDF1inds, EDF2inds, missingInds, removedOnsets1, removedOnsets2] = getEpochOnsets(pulseTimingPaths, systemTimeStart, teleporterEpochsPath, samplingRate)
% [epochOnsets, EDF1inds, EDF2inds] = getEpochOnsets(pulseTimingPaths, systemTimeStart, teleporterEpochsPath, samplingRate)
%
% Purpose: Take the timing of epochs onsets in Unix ticks and convert them
%   to indices of the EEG data.
%
% INPUT:
%   pulseTimingPaths: cell array of paths to the pulse timing files that
%       contain vectors of pulse times in Unix ticks and EEG indices
%   systemTimeStart: vector of epoch onset times in Unix ticks; output from
%       "sampleBehavioralData.m"
%   teleporterEpochsPath: path to the mat file containing the onsets of the
%       teleportation epochs in unix ticks, to be excluded from analysis;
%       in the case of one EDF, there is one variable (epochsEDF)
%   samplingRate: sampling rate of the EEG data in Hz
%
% OUTPUT:
%   epochOnsets: cell array of onset times in EEG indices, with onsets from
%       each EDF file in separate vectors
%   EDF1inds: indices of events that are contained within EDF1 (i.e., which
%       events from the regressors are in this EDF)
%   EDF2inds: same as above, but for EDF2
%   missingInds: indices of events that were lost in the break between EDF1
%       and EDF2
%   removedOnsets1: indices of the onsets in EDF1 that were removed because
%       they overlap with teleportation
%   removedOnsets2: indices of the onsets in EDF2 that were removed because
%       they overlap with teleportation
%
% Author: Lindsay Vass
% Date: 6 August 2015

warning off MATLAB:polyfit:RepeatedPointsOrRescale

% load pulse timing info
eegTime   = cell(length(pulseTimingPaths), 1);
unityTime = eegTime;
for thisEDF = 1:length(pulseTimingPaths)
    load(pulseTimingPaths{thisEDF});
    eegTime(thisEDF)   = {indEEG};
    unityTime(thisEDF) = {unityTicks};
end

% below is hard coded for 2 EDFs
fitRange = 5;
epochOnsets1 = [];
epochOnsets2 = [];
EDF1inds = [];
EDF2inds = [];
missingInds = [];
ind = 1;
for thisInterval = 1:length(systemTimeStart)
    thisStart = systemTimeStart(thisInterval);
    
    % trial before pulses started
    if thisStart < unityTime{1}(1)
        continue
        % trial in EDF1
    elseif thisStart < unityTime{1}(end)
        thisUnityTime = unityTime{1};
        thisEEGTime   = eegTime{1};
        epochOnsets1(end+1) = getEEGTime(thisStart, thisUnityTime, thisEEGTime, fitRange);
        EDF1inds(end + 1) = thisInterval;
        % trial lost between EDFs
    elseif thisStart < unityTime{2}(1)
        missingInds(end + 1) = thisInterval;
        continue
        % trial in EDF2
    else
        thisUnityTime = unityTime{2};
        thisEEGTime   = eegTime{2};
        epochOnsets2(end+1) = getEEGTime(thisStart, thisUnityTime, thisEEGTime, fitRange);
        EDF2inds(end + 1) = thisInterval;
    end
    
end

% load teleporter epochs info so we can exclude it later
load(teleporterEpochsPath);
teleEpochLength = table([1; 2], [1830; 2830], 'VariableNames', {'Index', 'LengthMs'});
epoch1Time      = table(eTime(1:min(missingEpochInds) - 1), epochsEDF1', 'VariableNames', {'Index', 'OnsetBins'});
epoch2Time      = table(eTime(max(missingEpochInds) + 1:length(eTime)), epochsEDF2', 'VariableNames', {'Index', 'OnsetBins'}); 
teleEpoch1Info  = innerjoin(teleEpochLength, epoch1Time);
teleEpoch2Info  = innerjoin(teleEpochLength, epoch2Time);

OffsetBins1     = table(round(teleEpoch1Info.OnsetBins + (teleEpoch1Info.LengthMs * (1/1000) * samplingRate)), 'VariableNames', {'OffsetBins'});
OffsetBins2     = table(round(teleEpoch2Info.OnsetBins + (teleEpoch2Info.LengthMs * (1/1000) * samplingRate)), 'VariableNames', {'OffsetBins'});
teleEpoch1Info  = [teleEpoch1Info, OffsetBins1];
teleEpoch2Info  = [teleEpoch2Info, OffsetBins2];

removedOnsets1 = [];
for thisTele = 1:length(teleEpoch1Info.OnsetBins)
    thisOnset  = teleEpoch1Info.OnsetBins(thisTele);
    thisOffset = teleEpoch1Info.OffsetBins(thisTele);
    firstOnset = find(thisOnset > epochOnsets1, 1, 'last' );
    lastOnset  = find(thisOffset < epochOnsets1, 1 ) - 1;
    removedOnsets1 = cat(1, removedOnsets1, [firstOnset:1:lastOnset]');
end

epochOnsets1(removedOnsets1) = [];



% determine the end time of the last missing epoch; if any of it is found
% at the beginning of EDF2, we'll remove it
missingOnset   = max(missingEpochs);
firstEDF2Onset = ticksEDF2(1);
onsetDiffBins  = round(((firstEDF2Onset - missingOnset) / 10000) * (1/1000) * samplingRate);
missingEpochLengthMs = teleEpochLength.LengthMs(find(teleEpochLength.Index == eTime(max(missingEpochInds))));
missingEpochLengthBins = round(missingEpochLengthMs * (1/1000) * samplingRate);
earliestValidOnset = epochsEDF2(1) - (onsetDiffBins - missingEpochLengthBins);

removedOnsets2 = [];
badOnsets = find(epochOnsets2 < earliestValidOnset);
removedOnsets2 = cat(1, removedOnsets2, badOnsets');
for thisTele = 1:length(teleEpoch2Info.OnsetBins)
    thisOnset  = teleEpoch2Info.OnsetBins(thisTele);
    thisOffset = teleEpoch2Info.OffsetBins(thisTele);
    firstOnset = find(thisOnset > epochOnsets2, 1, 'last' );
    lastOnset  = find(thisOffset < epochOnsets2, 1 ) - 1;
    removedOnsets2 = cat(1, removedOnsets2, [firstOnset:1:lastOnset]');
end

epochOnsets2(removedOnsets2) = [];

epochOnsets = cell(2,1);
epochOnsets{1} = epochOnsets1;
epochOnsets{2} = epochOnsets2;

end

function eegTime = getEEGTime(thisStart, thisUnityTime, thisEEGTime, fitRange)

% find pulses nearby our start time
timeDiff = abs(thisStart - thisUnityTime);
minTimeDiff = find(timeDiff == min(timeDiff));

% fit a linear model to get the time in EEG bins
fitInds = [minTimeDiff - fitRange:1:minTimeDiff + fitRange];
fitInds(fitInds < 1) = [];
fitInds(fitInds > length(thisUnityTime)) = [];
fit_P   = polyfit(thisUnityTime(fitInds), thisEEGTime(fitInds), 1);
eegTime = round(thisStart * fit_P(1) + fit_P(2));

end