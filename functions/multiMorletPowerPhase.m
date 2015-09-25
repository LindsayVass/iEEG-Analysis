function [frequencies, pow, pha] = multiMorletPowerPhase(eegData, frequencies, samplingRate, waveletCycles)
% multiMorletPowerPhase: return the power and phase for each frequency for
% each timepoint of EEG data after convolution with Morlet wavelets
% >> [frequencies, pow, pha] = multiMorletPowerPhase(eegData, frequencies, samplingRate, waveletCycles)
%
% Inputs:
%   eegData: vector of raw EEG data to be analyzed
%   frequencies: vector of frequencies in Hz at which to extract power and phase
%   samplingRate: samplingRate of eegData in Hz
%   waveletCycles: number of cycles of Morlet wavelets to use for
%       convolution (6 is recommended)
%
% Outputs:
%   frequencies: returns the same vector of frequencies
%   pow: frequencies x time points vector containing the power at each time
%       point for each frequency
%   pha: frequencies x time points vector containing the instantaneous
%       phase at each time point for each frequency
%
% Lindsay Vass
% 25 September 2015

pow = nan(length(frequencies), length(eegData));
pha = pow;
for thisFreq = 1:length(frequencies)
    [thisPow, thisPha] = morletPowerPhase(eegData, frequencies(thisFreq), samplingRate, waveletCycles);
    pow(thisFreq, :) = thisPow;
    pha(thisFreq, :) = thisPha;
end