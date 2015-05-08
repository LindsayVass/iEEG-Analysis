function [powerDist, freqs] = calcPowerDistribution(eegdata,samplerate,freqs,waveletwidth)
% CALCPOWERDISTRIBUTION - Calculate power distribution on entire EEG file.
%
% This function goes through the vector of EEG data for a single
% channel and calculates the power at each time point. It then fits a
% chi-square distribution to the wavelet power spectrum. Output is the
% distribution of power values at each frequency.
%
% FUNCTION: 
% calcPowerDistribution(eegdata,samplerate,freqs,waveletwidth)
%
% INPUT ARGS: 
%   eegdata = vector of EEG data
%   samplerate = sample rate in Hz of the EEG data
%   freqs = log-spaced set of frequencies at which power will be calculated
%   waveletwidth (optional) = number of wavelets for power calculation;
%   default is 6

if nargin <4
  waveletwidth = 6;
  if nargin <3
    error('eegdata, samplerate, and freqs are required to run');
  end
end

%we use some extra time at the beginning and end in order to avoid artifacts
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

fprintf('Calculating power...')
B=single(multienergyvec(eegdata,freqs,samplerate,waveletwidth));

% calc the mean fit to the background EEG spectrum
Blog = log10(double(B));
Pm = mean(Blog,2);
    
% get the fit
fprintf('Calc. fit...');

% IMPORTANT: chi_squarefit assumes that frequencies are logarithmically spaced!
[powerDist,R2] = chi_squarefit(freqs,Pm);
powerDist = powerDist';
    
% set the threshold to be the 95% percentile
% thresh = powerDist(:,10*amplitudeThresh+1);
    

function [B,t,f]=multienergyvec(S,f,Fs,width)
% function [B,t,f]=multienergyvec(S,f,Fs,width)
% compute the wavelet transform of the input data at the
% frequencies specified by f
% INPUT Args:
% s : signal
% f: vector with frequencies
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
% OUTPUT Args:
% B: matrix of length(freqs) by length(S) with the wavelet
% coefficients 
% t: time vector corresponding to the signal
% f: vector with frequencies for which power was computed
B = zeros(length(f),length(S));

fprintf('frequency ');
for a=1:length(f)
  fprintf('%d ',a);
  B(a,:)=energyvec(f(a),S,Fs,width);
end
fprintf('\n');

t = (1:size(S,2))/Fs;  



function [thresh,R2] = chi_squarefit(freqs, pows)
% function [thresh,R2] = chi_squarefit(freqs, pows)
% calculate the fit for the pepisode
nbins=1000;

pv=polyfit(log10(freqs),pows',1); % linear regression
means=10.^(polyval(pv,log10(freqs))); % these are the fitted mean powers

R=corrcoef(polyval(pv,log10(freqs)),pows');
R2=R(1,2)^2;

% Now, compute the chi2cdf vals.
thresh=[freqs;chi2inv(0:(1/nbins):(1-(1/nbins)),2)'*(means/2)];
% divide the means by two because the mean of the chi-square distribution is equal to the computed mean divided by the degrees of freedom




