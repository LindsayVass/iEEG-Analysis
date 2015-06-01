function [eegPower] = extractPowerWithBoundaries(EEG, chanInd, boundaries, frequencies, waveletWidth, bufferCycles)

% function [eegPower] = extractPower(EEG, chanInd,
% boundaries, frequencies, waveletWidth, bufferCycles)
% Extract the power at each time point of an EEG file that contains
% boundaries. Power will be estimated separately for each segment of the
% EEG dataset, where a segment is defined as the data between two boundary
% events. Before calculating power, buffers will be added to the beginning
% and end of the data, at a length corresponding to the time for the
% specified number of bufferCycles at the lowest frequency.
%
% INPUTS:
%   EEG: EEG structure from EEGLAB
%   chanInd: index of the channel in EEG.data to use for power estimation
%   boundaries: vector of boundary event onsets (in EEG samples)
%   frequencies: vector of frequencies in Hz at which to estimate power
%   waveletWidth (optional): number of Morlet wavelet cycles to use for
%       power estimation (default = 6)
%   bufferCycles (optional): number of cycles to use as a buffer for power
%       estimation (default = 3)
%
% OUTPUTS:
%   eegPower: struct containing for each segment the power at each
%       frequency

if nargin < 6
    bufferCycles = 3;
end

if nargin < 5
    waveletWidth = 6;
end

if nargin < 4
    error('Required arguments not supplied.')
end

% Setup output struct
eegPower = struct;

% make sure frequencies is a column vector
if size(frequencies, 2) == 1
else
    frequencies = frequencies';
end

fprintf(['\n\n\nSegment (' num2str(length(boundaries) + 1) ') total: ']);

% Loop through clean segments of the EEG
for thisSegment = 1:length(boundaries) + 1
    
    eegPower(thisSegment).segment = thisSegment;
    
    % Let the user know the current status
    fprintf([num2str(thisSegment) ' ']);
    
    [powerVal] = extractSegmentPower(EEG, chanInd, thisSegment, boundaries, frequencies, waveletWidth, bufferCycles);

    % Update struct
    eegPower(thisSegment).segment = struct('frequency',frequencies,'power',powerVal);
    
    
end % thisSegment

fprintf('\n\n');

end