function [binaryMatrix, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistributions, frequencies, eegData, samplingRate, amplitudeThresh, durationThresh)
% function [binaryMatrix, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistributions, frequencies, eegData, samplingRate, amplitudeThresh, durationThresh)
% Given a distribution of power values, determine whether there are
% sustained oscillations present in a given timeseries of EEG data. 
% To be counted as an oscillation, the power must exceed the value set by 
% amplitudeThresh and must continue for a period of time that exceeds 
% durationThresh. Power is calculated using 6-cycle Morlet wavelets and
% 3-cycle mirrored buffers are added to the eegData in order to estimate
% power near the data boundaries.
%
% INPUTS:
%   powerDistributions: structure of power distributions output from
%       calcPowerDistributionLKV.m
%   frequencies: vector of frequencies at which to assess pepisode
%   eegData: vector of EEG data
%   samplingRate: sampling rate of EEG data in Hz
%   amplitudeThresh (optional): amplitude threshold in percent of the
%       distribution (default = 95)
%   durationThresh (optional): duration threshold in cycles (default = 3)
%
% OUTPUTS:
%   binaryVector: logical matrix (frequencies x timepoints) that indicates 
%       for each time point whether there was a sustained oscillation 
%       present (1) or not (0)
%   percentTimePepisode: vector that indicates for each frequency the 
%       percent of time sustained oscillations were present in the supplied 
%       EEG data
%

if nargin < 6
    durationThresh = 3;
end

if nargin < 5
    amplitudeThresh = 95;
end

if nargin < 4
    error('Required arguments not supplied.')
end

% List the frequencies that we have power distributions for
freqList = cell2mat({powerDistributions.frequency});

% Check that the frequencies we're querying are the same frequencies we
% have distributions for
if sum(ismember(frequencies, freqList)) ~= length(frequencies)
    error('The frequencies you requested for pepisode calculations are not the frequencies provided in powerDistributions.')
    fprintf(['\nMissing frequencies: ' num2str(setdiff(frequencies, freqList)) '\n']);
end

% Buffer the EEG data
[bufferedData, bufferBins] = addMirroredBuffers(eegData, frequencies, samplingRate, 3);

% Extract power
powerVal = single(multienergyvec(bufferedData, frequencies, samplingRate, 6)); 

% Trim the buffers
powerVal = powerVal(:, bufferBins + 1:length(powerVal) - bufferBins);

% Initialize output matrix
binaryMatrix = zeros(length(frequencies), length(powerVal));

for thisFreq = 1:length(frequencies)
    
    freqInd = find(freqList == frequencies(thisFreq));
    powerDist = powerDistributions(freqInd).distribution;
    
    % Get the power threshold for this frequency
    powerThresh = powerDist(ceil(length(powerDist) * amplitudeThresh / 100));
    
    % Identify oscillations
    binaryVector = findPepisodes(powerVal(thisFreq,:), powerThresh, durationThresh, frequencies(thisFreq), samplingRate);
    
    % Add to the output matrix
    binaryMatrix(thisFreq, :) = binaryVector;
    
end % thisFreq

% Calculate percent pepisode
percentTimePepisode = mean(binaryMatrix, 2);

end