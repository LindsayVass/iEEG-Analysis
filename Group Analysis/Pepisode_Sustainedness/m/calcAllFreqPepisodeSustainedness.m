function [sustainednessTable] = calcAllFreqPepisodeSustainedness(powerDistributions, frequencies, eegData, samplingRate, eegTimes, amplitudeThresh, durationThresh)


% USAGE: [onsetMatrix, offsetMatrix, durationMatrix] = calcAllFreqPepisodeSustainedness(powerDistributions, frequencies, eegData, samplingRate, eegTimes, amplitudeThresh, durationThresh)
%
% Given a distribution of power values, determine whether there are
% sustained oscillations present in a given timeseries of EEG data.
% To be counted as an oscillation, the power must exceed the value set by
% amplitudeThresh and must continue for a period of time that exceeds
% durationThresh. Power is calculated using 6-cycle Morlet wavelets and
% 3-cycle mirrored buffers are added to the eegData in order to estimate
% power near the data boundaries. Will return three matrices indicating
% when oscillations started, when they stopped, and the duration.
%
% INPUTS:
%   powerDistributions: structure of power distributions output from
%       calcPowerDistributionLKV.m
%   frequencies: vector of frequencies at which to assess pepisode
%   eegData: vector of EEG data
%   samplingRate: sampling rate of EEG data in Hz
%   eegTimes: times (in seconds) of each bin of the epoch (can be found in
%       EEG.times of input data structure)
%   amplitudeThresh (optional): amplitude threshold in percent of the
%       distribution (default = 95)
%   durationThresh (optional): duration threshold in cycles (default = 2)
%
% OUTPUTS:
%   episodeMatrix: logical matrix (frequencies x timepoints) that indicates
%       for each time point whether there was a sustained oscillation
%       present (1) or not (0)
%   onsetMatrix: logical matrix (frequencies x timepoints) that indicates
%       for each time point whether an oscillation was initiated (1) or not
%       (0)
%   offsetMatrix: logical matrix (frequencies x timepoints) that indicates
%       for each time point whether an oscillation stopped (1) or not (0)

if nargin < 7
    durationThresh = 2;
end

if nargin < 6
    amplitudeThresh = 95;
end

if nargin < 5
    error('Required arguments not supplied.')
end

% Prepare the output table
sustainednessTable = table(nan(1), nan(1), nan(1), nan(1));

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

for thisFreq = 1:length(frequencies)
    
    freqInd = find(freqList == frequencies(thisFreq));
    powerDist = powerDistributions(freqInd).distribution;
    
    % Get the power threshold for this frequency
    powerThresh = powerDist(ceil(length(powerDist) * amplitudeThresh / 100));
    
    % Identify oscillations
    [onsetTime, offsetTime, duration] = calcPepisodeSustainedness(powerVal(thisFreq, :), powerThresh, durationThresh, frequencies(thisFreq), samplingRate, eegTimes);
    
    % Add to the output table
    if isempty(onsetTime)
        continue
    else
        tempTable = table(onsetTime', offsetTime', duration', repmat(frequencies(thisFreq), length(onsetTime), 1));
        sustainednessTable = [sustainednessTable; tempTable];
    end % thisFreq
    
end

sustainednessTable.Properties.VariableNames = {'Onset', 'Offset', 'Duration', 'Frequency'};
sustainednessTable(1,:) = [];