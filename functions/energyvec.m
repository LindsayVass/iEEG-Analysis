function powerVals = energyvec(frequency,eegData,samplingRate,waveletWidth)
% function y = energyvec(f,s,Fs,width)
%
% Return a vector containing the power as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% 
% INPUTS:
%   frequency: frequency in Hz at which to estimate power
%   eegData: vector of EEG data
%   samplingRate: sampling rate in Hz
%   waveletWidth : width of Morlet wavelet in cycles (>= 5 suggested).
%
% 

dt = 1/samplingRate;
sf = frequency/waveletWidth;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(frequency,t,waveletWidth);

powerVals = conv(eegData,m);

powerVals = abs(powerVals).^2;
powerVals = powerVals(ceil(length(m)/2):length(powerVals)-floor(length(m)/2));