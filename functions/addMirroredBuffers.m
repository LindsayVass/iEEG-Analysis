function [bufferedData, bufferBins] = addMirroredBuffers(eegData, frequencies, samplingRate, bufferCycles)
% function [bufferedData, bufferBins] = addMirroredBuffers(eegData, frequencies, samplingRate, bufferCycles)
% Add mirrored version of EEG signal to beginning and end of signal as a
% buffer in order to calculate power at time points near the edge of the
% signal. The length of the buffer is equivalent to the time it takes to go
% through bufferCycles at the lowest frequency. For example, for
% min(frequencies) = 1 and bufferCycles = 3, the buffer is 3000 ms.
%
% INPUTS:
%   eegData: vector of EEG data
%   frequencies: vector of frequencies in Hz
%   samlpingRate: sampling rate in Hz
%   bufferCycles (optional): number of cycles to use for buffer length
%       (default = 3)
%
% OUTPUTS:
%   bufferedData: vector of EEG data with buffers added to beginning and
%       end
%   bufferBins: length of the buffer in EEG samples; returned so you can
%       easily remove the buffer later

if nargin < 4
    bufferCycles = 3;
end

if nargin < 3
    error('Required arguments not supplied.')
end

% Calculate the lengthe of the buffer in milliseconds and in EEG samples
bufferMS   = bufferCycles * 1 / min(frequencies) * 1000;
bufferBins = round(bufferMS * 1/1000 * samplingRate);
numBuffers = ceil(bufferBins / length(eegData));

bufferData = [];
for thisChunk = 1:numBuffers
    
    % Mirror the signal on every other chunk so it's continuous
    if mod(thisChunk, 2) == 1
        bufferData = cat(2, bufferData, fliplr(eegData));
    else
        bufferData = cat(2, bufferData, eegData);
    end
    
end % thisChunk

% Trim excess data from buffer
bufferData = bufferData(1:bufferBins);

% Return buffered data
bufferedData = [fliplr(bufferData), eegData, bufferData];


end