function [power,times,frequencies] = multienergyvec(eegData,frequencies,samplingRate,waveletWidth)
% function [B,t,f] = multienergyvec(eegData,f,Fs,width)
% compute the wavelet transform of the input data at the
% frequencies specified by f
%
% INPUTS:
%   eegData: signal
%   frequecies: vector of frequencies in Hz at which to perform wavelet 
%       convolution
%   samplingRate: sampling rate in Hz
%   waveletWidth : width of Morlet wavelet in cycles (>= 5 suggested).
% 
% OUTPUTS:
%   power: matrix of length(frequencies) by length(eegData) with the wavelet
%       coefficients 
%   times: time vector corresponding to the signal time in seconds
%   frequencies: vector with frequencies for which power was computed

power = zeros(length(frequencies),length(eegData));


for thisFrequency = 1:length(frequencies)
  power(thisFrequency,:)=energyvec(frequencies(thisFrequency),eegData,samplingRate,waveletWidth);
end

times = (1:size(eegData,2))/samplingRate; 