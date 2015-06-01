function [powerVal] = extractSegmentPower(EEG, chanInd, segmentNumber, boundaries, frequencies, waveletWidth, bufferCycles)
% function [powerVal] = extractSegmentPower(EEG, chanInd, segmentNumber,
% boundaries, frequencies, waveletWidth, bufferCycles)
% Extract the power at each time point within a given segment of the EEG,
% where a segment is the data between two boundaries. Before calculating
% power, buffers will be added to the beginning and end of the data, at a
% length corresponding to the time for the specificed number of
% bufferCycles at the lowest frequency.
%
% INPUTS:
%   EEG: EEG structure from EEGLAB
%   chanInd: index of the channel in EEG.data to use for power estimation
%   segmentNumber: which segment of the EEG data to extract power from
%   boundaries: vector of boundary onsets in EEG samples
%   frequencies: vector of frequencies in Hz at which to estimate power
%   waveletWidth (optional): number of Morlet wavelet cycles to use for
%       power estimation (default = 6)
%   bufferCycles (optional): number of cycles to use as a buffer for power
%       estimation (default = 3)
%
% OUTPUT:
%   powerVal: matrix (frequencies x timepoints) of power values at each
%       frequency and timepoint

if nargin < 7
    bufferCycles = 3;
end

if nargin < 6
    waveletWidth = 6;
end

if nargin < 5
    error('Required arguments not supplied.')
end

powerVal = [];

% Specify the first bin of the segment. If this if the first segment, then
% the beginning of the EEG file is the first bin.
if segmentNumber == 1
    
    firstBin = 1;
    
    % If there was a problem at the beginning of the recording,
    % (usually beginning of 2nd EEG), the first boundary will
    % be zero so skip it.
    if boundaries(segmentNumber) == 0
        return
    end
    
else
    firstBin = boundaries(segmentNumber - 1) + 1;
end

% Specify the last bin of the segment. If this is the last segment,
% the last bin is the length of the EEG.
if segmentNumber == length(boundaries) + 1
    lastBin = size(EEG.data, 2);
else
    lastBin = boundaries(segmentNumber);
end

% Extract EEG data for this segment for this channel
eegData = EEG.data(chanInd, firstBin:lastBin);

% Add buffers to EEG data
[bufferedData, bufferBins] = addMirroredBuffers(eegData, frequencies, EEG.srate, bufferCycles);

% Extract power
powerVal = single(multienergyvec(bufferedData, frequencies, EEG.srate, waveletWidth)); 

% Remove buffers from data
powerVal = powerVal(:, bufferBins + 1:length(powerVal) - bufferBins);

end