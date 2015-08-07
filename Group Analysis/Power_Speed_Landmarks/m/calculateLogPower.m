function [logPowerVals, EEG] = calculateLogPower(eegPath, chanName, frequencies, numWavelets, bufferCycles)
% [logPowerVals] = calculateLogPower(eegPath, chanName, frequencies, numWavelets, bufferCycles)
%
% Purpose: Calculate the log(Power) at a set of frequencies for a given
%   electrode. The data at the beginning and end will be mirrored to
%   produce a buffer that allows us to accurately estimate power at the
%   edge of the recordings. The duration, set by bufferCycles (optional
%   input) defaults to 3 cycles at the lowest frequency.
%
% INPUT:
%   eegPath: path to the unepoched EEG data
%   chanName: string containing the name of the desired electrode
%   frequencies: vector of frequencies in Hz at which to calculate power
%
% OPTIONAL INPUT:
%   numWavelets: number of wavelets to use for convolution to estimate
%       power (default = 6)
%   bufferCycles: number of cycles to use for the buffer at the edge of the
%       recordings (default = 3)
%
% OUTPUT:
%   logPowerVals: frequency x timepoint matrix of log(Power) values
%   EEG: the EEGLAB structure on which the analysis was run
%
% Author: Lindsay Vass
% Date: 6 August 2015

if nargin < 5
    bufferCycles = 3;
    if nargin < 4
        numWavelets = 6;
    end
end

EEG = pop_loadset(eegPath);

% find the index of this electrode
chanInd = find(strcmpi(chanName, {EEG.chanlocs.labels}));

% initialize output
logPowerVals = nan(length(frequencies), size(EEG.data, 2));

% extract data for this electrode
rawData = EEG.data(chanInd,:);

% add buffers
[bufferedData, bufferBins] = addMirroredBuffers(rawData, frequencies, EEG.srate, bufferCycles);

% calculate power
[powerVals, ~, ~] = multienergyvec(bufferedData, frequencies, EEG.srate, numWavelets);

% trim the data to remove the buffers
powerValsTrimmed = powerVals(:, bufferBins + 1:length(powerVals) - bufferBins);

% take the log
logPowerVals = log(powerValsTrimmed);

end