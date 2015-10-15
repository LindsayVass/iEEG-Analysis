function [binaryMatrix, meanPercentTimePepisode, meanPower] = calcPepisodeWholeExpt(powerDistributions, frequencies, eegList, chanName, amplitudeThresh, durationThresh)
% >> [binaryMatrix, percentTimePepisode] = calcPepisodeWholeExpt(powerDistributions, frequencies, eegList, chanName, amplitudeThresh, durationThresh)
%
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
%   eegList: cell array of paths to the unepoched EEG data
%   chanName: char of the electrode name
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

allBinaryMatrices = cell(1, length(eegList));
allPercentTimePepisode = cell(1, length(eegList));
allPowerVal = cell(1, length(eegList));

% Loop through EEG datasets
for thisEEG = 1:length(eegList)
    
    % Load the EEG data
    EEG = pop_loadset(eegList{thisEEG});
    
    % Identify the boundaries within the EEG file
    boundaries = findBoundaries(EEG);
    
    % find the index of this channel
    chanInd = find(strcmpi(chanName, {EEG.chanlocs.labels}));
    
    % Calculate power
    [powerVal] = extractPowerWithBoundaries(EEG, chanInd, boundaries, frequencies, 6);
    
    % Concatenate into one matrix and take the mean
    tempPowerVal = cell(1, length(powerVal));
    for i = 1:length(powerVal)
       tempPowerVal{i} = powerVal(i).segment.power;
    end
    tempPowerVal = [tempPowerVal{:}];
    meanPowerVal = mean(tempPowerVal, 2);
    
    % Initialize output matrix
    binaryMatrix = cell(length(frequencies), length(powerVal));
    
    for thisSegment = 1:length(powerVal)
        
        thisData = powerVal(thisSegment).segment.power;
        
        if (size(thisData, 1) == 0)
            continue
        end
        
        for thisFreq = 1:length(frequencies)
            
            freqInd = find(freqList == frequencies(thisFreq));
            powerDist = powerDistributions(freqInd).distribution;
            
            % Get the power threshold for this frequency
            powerThresh = powerDist(ceil(length(powerDist) * amplitudeThresh / 100));
            
            % Identify oscillations
            binaryVector = findPepisodes(thisData(thisFreq,:), powerThresh, durationThresh, frequencies(thisFreq), EEG.srate);
            
            % Add to the output matrix
            binaryMatrix{thisFreq, thisSegment} = binaryVector;
            
        end % thisFreq
        
    end % thisSegment
    
    % Concatenate all the binary data
    concatBinaryMatrix = [];
    for i = 1:size(binaryMatrix, 1)
       concatBinaryMatrix(i, :) = [binaryMatrix{i,:}]; 
    end
    
    
    allBinaryMatrices{thisEEG} = concatBinaryMatrix;
    allPowerVal{thisEEG} = tempPowerVal;
    
    % Calculate percent pepisode
    percentTimePepisode = mean(concatBinaryMatrix, 2);
    
    allPercentTimePepisode{thisEEG} = percentTimePepisode;
    
end % thisEEG

% Take average across all EEGs
concatPepisode = [allPercentTimePepisode{:}];
meanPercentTimePepisode = mean(concatPepisode, 2);

concatPower = [allPowerVal{:}];
meanPower = mean(concatPower, 2);

end