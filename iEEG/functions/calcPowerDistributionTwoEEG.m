function [powerDist, freqs] = calcPowerDistributionTwoEEG(eegdata1,eegdata2,samplerate,freqs,waveletwidth)
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

if nargin <5
  waveletwidth = 6;
  if nargin <4
    error('eegdata1, eegdata2, samplerate, and freqs are required to run. If you only have one EEG file, use "calcPowerDistribution" instead.');
  end
end

%we use some extra time at the beginning and end in order to avoid artifacts
shoulderMS = 500;
shoulder = round(shoulderMS*samplerate/1000); % shoulder in samples

fprintf('Calculating power...')
B1=single(multienergyvec(eegdata1,freqs,samplerate,waveletwidth));
B2=single(multienergyvec(eegdata2,freqs,samplerate,waveletwidth));
B = cat(2, B1, B2);
    
% calc the mean fit to the background EEG spectrum
Blog = log10(double(B));
Pm = mean(Blog,2);
    
% get the fit
fprintf('Calc. fit...');
% IMPORTANT: chi_squarefit assumes that frequencies are logarithmically spaced!
[powerDist,R2] = chi_squarefit(freqs,Pm);
powerDist = powerDist';
    
% set the threshold to be the 95% percentile
% thresh = all(:,10*amplitudeThresh+1);
 
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




% [P,holes,detected]=episodeid(timecourse,powerthresh,durthresh,shoulder,plotit)
%
% This script identifies oscillatory episodes where the threshold is given
% leadno can be a vector of leads to analyse.
%
% Parameters:
% timecourse - the wavelet transformed signal to analyse
% powerthresh - wavelet power threshold
% durthresh- minimum time (duration threshold)
% shoulder- optionally, you can include a shoulder so that episodes
%           identification won't be subject to edge artifacts. This value
%           should be specified in _samples_
% [plotit]- optionally set this to 1 to have the script show you how it's
%           working.
%              Red - the wavelet power timecourse you passed in
%                    'timecourse'
%              Green - the 'detected' variable (times where episodes were
%                      found)
%              Cyan - the shoulder (which are excluded from the Pepisode
%                     calculation, but included while detecting episodes)
%
% Returns:
% P- Pepisode (percentage of time occupied by oscillatory episodes)
% holes- start and end times of detected episodes. It has units of
%        seconds.
% detected- a binary vector containing a 1 where an episode was detected
%           and a 0 otherwise

function detected=episodeid(timecourse,powerthresh,durthresh,shoulder,plotit)


TIMECOURSE=timecourse;
if(nargin<5) plotit=0; end
NAN_CODE=-99; samplerate=eegparams('samplerate');
t=(1:length(timecourse))/samplerate;
THRESH=powerthresh*ones(size(timecourse));
timecourse(find(timecourse<THRESH))=0; % zero anything under the threshold, same threshold

detected=zeros(1,length(timecourse)); % to store all episodes

x=(timecourse>0); % a binary vector: pass or not

dx=diff(x); pos=find(dx==1)+1; neg=find(dx==-1)+1; % show the +1 and -1 edges
clear dx;
% now do all the special cases to handle the edges
if(isempty(pos) & isempty(neg))
  if(find(x)>0) H=[1;length(timecourse)]; else H=[]; end % all episode or none
elseif(isempty(pos)) H=[1;neg]; % i.e., starts on an episode, then stops
elseif(isempty(neg)) H=[pos;length(timecourse)]; % starts, then ends on an ep.
else
  if(pos(1)>neg(1))
      pos=[1 pos];
  end; % we start with an episode
  if(neg(length(neg))<pos(length(pos)))
      neg=[neg length(timecourse)];
  end; % we end with an episode
  H=[pos;neg]; % NOTE: by this time, length(pos)==length(neg), necessarily
end; % special-casing, making the H double-vector
clear x pos neg;
if(~isempty(H)) % more than one "hole"
  % find epochs lasting longer than minNcycles*period
  goodep=find((H(2,:)-H(1,:))>=durthresh);
  if(isempty(goodep))
  	H=[];
  else
  	H=H(:,goodep);
	% fprintf(1,'ge ');
  end;
  % this becomes detected vector
  for h=1:size(H,2) detected(H(1,h):H(2,h))=1; end;
end % more than one "hole"

if(shoulder~=0) % handle shoulder
  detected(1:shoulder)=0;
  detected((length(detected)-shoulder):length(detected))=0;
end

% now, consolidate all the intervals

x=(detected>0); % a binary vecor: pass or not
dx=diff(x); pos=find(dx==1)+1; neg=find(dx==-1)+1; % show the +1 and -1 edges
clear dx;
% now do all the special cases to handle the edges
if(isempty(pos) & isempty(neg))
  if(find(x)>0) holes=[1;length(timecourse)]; else holes=[]; end % all episode or none
elseif(isempty(pos)) holes=[1;neg]; % i.e., starts on an episode, then stops
elseif(isempty(neg)) holes=[pos;length(timecourse)]; % starts, then ends on an ep.
else
  if(pos(1)>neg(1))
  	pos=[1 pos]; end; % we start with an episode
  % we end with an episode
  if(neg(length(neg))<pos(length(pos)))
  	neg=[neg length(timecourse)]; end;
  holes=[pos;neg]; % NOTE: by now, length(pos)==length(neg), necessarily
end; % special-casing, making the H double-vector
clear x pos neg;

if(plotit)
  plot((1:length(TIMECOURSE))/samplerate,TIMECOURSE,'r-');
  axis('tight'); replot; ax=axis; hold on;
  miny=min(TIMECOURSE); maxy=max(TIMECOURSE); yval=miny+0.2*(maxy-miny);
  plot(holes/samplerate,yval*ones(size(holes)),'+-g');
  legend('off'); xlabel('Time [s]'); ylabel('Wavelet Power(t)');
  yval=miny+0.1*(maxy-miny);
  plot([1 shoulder]/samplerate,[yval yval],'co-');
  plot((length(TIMECOURSE)-[0 shoulder])/samplerate,[yval yval],'co-');
  hold off;
end % plotit

u=find(detected);
if(isempty(u)), P=0; else, P=length(u)/(length(detected)-2*shoulder); end
