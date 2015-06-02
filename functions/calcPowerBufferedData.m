function powerVal = calcPowerBufferedData(bufferedData, bufferBins, frequencies, samplingRate, waveletWidth)

% function powerVal = calcPowerBufferedData(bufferedData, bufferBins, frequencies, samplingRate, waveletWidth)
% Extract power at each timepoint of data using Morlet wavelets and then
% trim the buffers.
%
% INPUTS: 
%   bufferedData: vector of EEG data with buffers added to beginning and
%       end (see addMirroredBuffers.m)
%   bufferBins: size of the buffer period in matrix elements
%   frequencies: vector of frequencies in Hz at which to extract power
%   samplingRate: sampling rate of the bufferedData in Hz
%   waveletWidth (optional): number of cycles to use for Morlet wavelet
%       convolution (default = 6)
%
% OUTPUT:
%   powerVal: matrix (frequencies x timepoints) of power values at each
%       frequency and timepoint
%

if nargin < 5
    waveletWidth = 6;
end

if nargin < 4
    error('Required arguments not supplied.')
end

% Extract power
powerVal = single(multienergyvec(bufferedData, frequencies, samplingRate, waveletWidth)); 

% Remove buffers from data
powerVal = powerVal(:, bufferBins + 1:length(powerVal) - bufferBins);

end