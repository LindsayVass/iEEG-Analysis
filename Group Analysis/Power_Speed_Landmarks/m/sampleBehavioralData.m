function [systemTimeStart, systemTimeEnd, landmarkTypeSamp, speedSamp] = sampleBehavioralData(systemTime, xPos, zPos, intervalMs)
% [systemTimeStart, systemTimeEnd, landmarkTypeSamp, speedSamp] = sampleBehavioralData(systemTime, xPos, zPos, intervalMs)
% 
% Purpose: Take the output of "parseBehavioralTxt.m" and divide it into
%   intervals whose length in milliseconds is defined by intervalMs. Then,
%   return the mean (for numeric) or modal (for string) value for each
%   interval.
%
% INPUT:
%   systemTime: vector of time values in Unix ticks
%   xPos: vector of x coordinates
%   zPos: vector of z coordinates
%   intervalMs: length of the desired interval in milliseconds
%
% OUTPUT:
%   systemTimeStart: vector of values representing the onset of each
%       interval in Unix ticks
%   systemTimeEnd: vector of values representing the offset of each
%       interval in Unix ticks
%   landmarkTypeSamp: cell array containing the modal landmark type across
%       the interval; landmark types include "RICH", "POOR", and "CENTRAL",
%       which are defined based on the avatar's xz coordinate
%   speedSamp: vector of values representing the avatar's mean speed across
%       the interval
%
% Author: Lindsay Vass
% Date: 5 August 2015

% get sampling rate of the text file
txtIntervalMs      = round(mean(diff(systemTime))) / 10000; % 10000 ticks per millisecond
samplesPerInterval = round(intervalMs / txtIntervalMs);

% make vector to split the data
splitInd = repmat([1:1:length(systemTime)/samplesPerInterval], [samplesPerInterval, 1]);
splitInd = reshape(splitInd, [], 1);

% for each sample, calculate mean speed and modal landmark richness
speedSamp        = nan(max(splitInd), 1);
landmarkTypeSamp = cell(max(splitInd), 1);
systemTimeStart  = speedSamp;
systemTimeEnd    = speedSamp;

for i = 1:max(splitInd)
    
    % get start and end times
    timeSamp = systemTime(splitInd == i);
    systemTimeStart(i) = timeSamp(1);
    systemTimeEnd(i)   = timeSamp(end);
    
    % calculate speed for this interval
    xSamp = xPos(splitInd == i);
    zSamp = zPos(splitInd == i);
    speedSamp(i) = calcSpeed(xSamp, zSamp, txtIntervalMs);
    
    % determine the modal landmark richness
    landmarkTypes = getLandmarkType(xSamp, zSamp);
    landmarkTypeSamp{i} = getModalLandmark(landmarkTypes);
    
end

end

function meanSpeed = calcSpeed(xSamp, zSamp, txtIntervalMs)

distMat = dist([cat(1, xSamp', zSamp')]);

% restrict to distances between consecutive points
p1 = [1:1:length(xSamp)-1];
p2 = [2:1:length(xSamp)];
distances = distMat(sub2ind(size(distMat), p1, p2));
meanSpeed = mean(distances / txtIntervalMs);

end

function landmarkTypes = getLandmarkType(xSamp, zSamp)

landmarkTypes = cell(length(xSamp), 1);

for j = 1:length(xSamp)
    
    thisX = xSamp(j);
    thisZ = zSamp(j);
    
    % RICH arms
    if (thisZ > 586 | thisZ < 488)
        landmarkTypes(j) = {'RICH'};
    % POOR arms    
    elseif (thisX > 536 | thisX < 438)
        landmarkTypes(j) = {'POOR'};
    % CENTRAL plaza
    else
        landmarkTypes(j) = {'CENTRAL'};
    end
    
end

end

function modalLandmark = getModalLandmark(landmarkTypes)

stringList  = unique(landmarkTypes);
stringCount = zeros(length(stringList), 1);

for k = 1:length(stringList)
    stringCount(k) = length(find(strcmpi(stringList{k}, landmarkTypes)));
end

[~, ind] = max(stringCount);
modalLandmark = stringList(ind);

end