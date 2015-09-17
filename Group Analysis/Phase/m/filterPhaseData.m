function [phaseDataSample, goodFreq, badFreq, timeMs] = filterPhaseData(phaseData, timeInterval, eegTimes, frequencies, samplingRate, numCycles)
% filterPhaseData: Take a dataset of phases and filter data, keeping only
% values within the specified time interval, and at frequencies that exceed
% the minimum determined by numCycles (must have at least numCycles within
% the time interval)
% >> [phaseDataSample, goodFreq, badFreq] = filterPhaseData(phaseData, timeInterval, eegTimes, frequencies, samplingRate, numCycles)
%
% Inputs:
%   phaseData: electrodes x timepoints x trials x frequencies dataset of
%       phases
%   timeInterval: two-value vector indicating start and end times of time
%       interval in milliseconds
%   eegTimes: vector from EEG.times
%   frequencies: vector of frequencies in Hz
%   samplingRate: sampling rate of the EEG data in Hz
%   numCycles: number of cycles of data required for phase analysis; if a
%       frequency is too slow to achieve numCycles within timeInterval, it
%       will be excluded
%
% Outputs:
%   phaseDataSample: electrodes x timepoints x trials x frequencies dataset
%       of phases
%   goodFreq: list of the frequencies included in the phaseDataSample
%   badFreq: list of excluded frequencies
%   timeMs: vector of timepoints in ms (0 = teleporter entry)
%
% Lindsay Vass
% 15 September 2015

timeInds  = find(eegTimes >= timeInterval(1) & eegTimes <= timeInterval(2));
timeMs  = eegTimes(timeInds);
minFreq   = numCycles / (length(timeInds) / samplingRate);
goodFreq  = frequencies(frequencies >= minFreq);
badFreq   = frequencies(frequencies < minFreq);

phaseDataSample = phaseData(:, timeInds, :, frequencies >= minFreq);