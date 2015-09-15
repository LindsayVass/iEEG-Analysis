function [pow, pha] = morletPowerPhase(eegData, frequency, samplingRate, waveletCycles)
% morletPowerPhase: return the power and phase for each timepoint of EEG
% data after convolution with Morlet wavelets
% >> [pow, pha] = morletPowerPhase(eegData, frequency, samplingRate, waveletCycles)
%
% Inputs:
%   eegData: vector of raw EEG data to be analyzed
%   frequency: frequency in Hz at which to extract power and phase
%   samplingRate: samplingRate of eegData in Hz
%   waveletCycles: number of cycles of Morlet wavelets to use for
%       convolution (6 is recommended)
%
% Outputs:
%   pow: vector of the same size as eegData containing the power at each
%       time point for the selected frequency
%   pha: vector of the same size as eegData containing the instantaneous
%       phase at each time point for the selected frequency
%
% Lindsay Vass
% 15 September 2015

% FFT parameters (use next-power-of-2)
time = -1:1/samplingRate:1;
halfWaveletSize = (length(time)-1)/2;
nWavelet      = length(time);
nData         = length(eegData);
nConvolution  = nWavelet+nData-1;
nConvPow2    = pow2(nextpow2(nConvolution));

fftData = fft(eegData,nConvPow2);

wavelet    = (pi*frequency*sqrt(pi))^-.5 * exp(2*1i*pi*frequency.*time) .* exp(-time.^2./(2*( waveletCycles /(2*pi*frequency))^2))/frequency;
fftWavelet = fft(wavelet,nConvPow2);

convolutionResultFft = ifft(fftWavelet.*fftData,nConvPow2);
convolutionResultFft = convolutionResultFft(1:nConvolution); % note: here we remove the extra points from the power-of-2 FFT
convolutionResultFft = convolutionResultFft(halfWaveletSize+1:end-halfWaveletSize);

pow = abs(convolutionResultFft).^2;
pha = angle(convolutionResultFft);