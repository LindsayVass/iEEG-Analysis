function [epochOnsets, EDF1inds, EDF2inds, missingInds] = getEpochOnsets(pulseTimingPaths, systemTimeStart)
% [epochOnsets, EDF1inds, EDF2inds] = getEpochOnsets(pulseTimingPaths, systemTimeStart)
%
% Purpose: Take the timing of epochs onsets in Unix ticks and convert them
%   to indices of the EEG data.
%
% INPUT:
%   pulseTimingPaths: cell array of paths to the pulse timing files that
%       contain vectors of pulse times in Unix ticks and EEG indices
%   systemTimeStart: vector of epoch onset times in Unix ticks; output from
%       "sampleBehavioralData.m"
%
% OUTPUT:
%   epochOnsets: cell array of onset times in EEG indices, with onsets from
%       each EDF file in separate vectors
%   EDF1inds: indices of events that are contained within EDF1 (i.e., which
%       events from the regressors are in this EDF)
%   EDF2inds: same as above, but for EDF2
%   missingInds: indices of events that were lost in the break between EDF1
%       and EDF2
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
epochOnsets = cell(size(eegTime));
ind1 = 1;
ind2 = 1;
EDF1inds = [];
EDF2inds = [];
missingInds = [];
for thisInterval = 1:length(systemTimeStart)
    thisStart = systemTimeStart(thisInterval);
    
    % trial before pulses started
    if thisStart < unityTime{1}(1)
        continue
        % trial in EDF1
    elseif thisStart < unityTime{1}(end)
        thisUnityTime = unityTime{1};
        thisEEGTime   = eegTime{1};
        epochOnsets{1}(ind1) = getEEGTime(thisStart, thisUnityTime, thisEEGTime, fitRange);
        ind1 = ind1 + 1;
        EDF1inds(end + 1) = thisInterval;
        % trial lost between EDFs
    elseif thisStart < unityTime{2}(1)
        missingInds(end + 1) = thisInterval;
        continue
        % trial in EDF2
    else
        thisUnityTime = unityTime{2};
        thisEEGTime   = eegTime{2};
        epochOnsets{2}(ind2) = getEEGTime(thisStart, thisUnityTime, thisEEGTime, fitRange);
        ind2 = ind2 + 1;
        EDF2inds(end + 1) = thisInterval;
    end
    
    
end

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