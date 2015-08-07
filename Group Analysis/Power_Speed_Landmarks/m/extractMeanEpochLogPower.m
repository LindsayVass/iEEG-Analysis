function [meanLogPowerVals] = extractMeanEpochLogPower(logPowerVals, epochOnsets, epochLengthMs, samplingRate)
% [meanLogPowerVals] = extractMeanEpochLogPower(logPowerVals, epochOnsets, epochLengthMs)
%
% Purpose: For each epoch, calculate the mean log(Power) at each frequency.
%
% INPUT:
%   logPowerVals: frequencies x time points vector of log(Power) values
%   epochOnsets: vector of indices indicating the start time for each epoch
%   epochLengthMs: length of the epoch in milliseconds
%   samplingRate: sampling rate of the EEG data in Hz
%
% OUTPUT:
%   epochedLogPowerVals: frequency x epoch matrix of mean log(Power) values
%
% Author: Lindsay Vass
% Date: 6 August 2015

epochLengthBins = round(epochLengthMs * (1/1000) * samplingRate) - 1; % subtract 1 because epochOnset will be the first bin

meanLogPowerVals = nan(size(logPowerVals, 1), length(epochOnsets));
for thisEpoch = 1:length(epochOnsets)
    thisStart = epochOnsets(thisEpoch);
    thisEnd   = epochOnsets(thisEpoch) + epochLengthBins;
    thisData = logPowerVals(:, thisStart:thisEnd);
    meanLogPowerVals(:, thisEpoch) = mean(thisData, 2);
end

end
