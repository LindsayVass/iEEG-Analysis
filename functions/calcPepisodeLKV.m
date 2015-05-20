function calcPepisodeLKV(eegList, chanNames, saveStem, frequencies, durationThresh, amplitudeThresh)
% calcPepisodeLKV: Calculate pepisode on the clean segments of the EEG.
% This function will loop through the channels in "chanNames" and calculate
% the pepisode transform of the EEG data. The function will first loop
% through all EEG files and calculate the power at each time point. It will
% then concatenate all the power values and fit a chi-square distribution
% to the mean power across frequencies. The function will then create a
% binary vector that indicates for each time point whether the power at
% that time exceeded both the duration threshold (default = 3 cycles) and
% the amplitude threshold (default = 95% of the chi-square distribution). 
%
% INPUTS:
%   eegList: cell array of paths to the EEG(s) to analyze
%
%   chanNames: cell array of channel names to analyze
%
%   frequencies: vector of log-spaced frequencies to use for analysis
%
%   saveStem: string that includes the path to the folder to store the file
%       as well as the initial segment of the fileName (the script will
%       append the electrode name); for example
%       '/path/to/pepisode/UCDMC15_TeleporterA_pepisode_'
%
%   durationThresh (optional): duration threshold in number of cycles
%
%   amplitudeThresh (optional): amplitude threshold in percentage of power
%       distribution; default = 95

if nargin < 6
  amplitudeThresh = 95;
  if nargin < 5
    durationThresh = 3;
    if nargin < 4
        error('Not enough arguments supplied.')
    end
  end
end

% wavelet width
width = 6;

% initialize eeglab
eeglab;

% Loop over channels
for thisChan = 1:length(chanNames)
    
    % Let the user know the current status
    fprintf(['\n\n\nCalculating power for channel ' num2str(thisChan) ' of ' num2str(length(chanNames)) '\n']);
    
    % Initialize a cell array to hold all of the power values from each segment
    % of each EEG
    powerHolder = cell(size(eegList));
    
    % Initialize a cell array to hold the duration of each segment of each EEG
    durationHolder = cell(size(eegList));
    
    % Loop over EEGs
    for thisEEG = 1:length(eegList)
        
        % Initialize a vector to hold all power values from this EEG
        eegPower = [];
        
        % Initialize a vector to hold all segment durations from this EEG
        segmentDurations = [];
        
        % load EEG
        EEG = pop_loadset(eegList{thisEEG});
        
        % we use some extra time at the beginning and end in order to avoid artifacts
        shoulderMS = 500;
        shoulder = round(shoulderMS * EEG.srate / 1000); % shoulder in samples

        
        % find the index of this channel
        chanInd = find(strcmpi(chanNames{thisChan}, {EEG.chanlocs.labels}));
        
        % Get the indices of the boundaries in the file
        boundaries = floor(cell2mat({EEG.event.latency}));
        
        % Let the user know the current status
         fprintf(['EEG ' num2str(thisEEG) ' of ' num2str(length(eegList)) '\nSegment (' num2str(length(boundaries) + 1) ' total):\n']);
        
        
        % Loop through clean segments of the EEG
        for thisSegment = 1:length(boundaries) + 1
            
            % Let the user know the current status
            fprintf([num2str(thisSegment) ' ']);
            
            % Specify the first bin of the segment. If this is the first
            % segment, the first bin is 1. 
            if thisSegment == 1
                
                firstBin = 1;
                
                % If there was a problem at the beginning of the recording,
                % (usually beginning of 2nd EEG), the first boundary will
                % be zero so skip it.
                if boundaries(thisSegment) == 0
                    continue
                end
                
            else
                firstBin = boundaries(thisSegment - 1) + 1;
            end
            
            % Specify the last bin of the segment. If this is the last segment,
            % the last bin is the length of the EEG. 
            if thisSegment == length(boundaries) + 1
                lastBin = size(EEG.data, 2);
            else
                lastBin = boundaries(thisSegment);
            end
            
            % Save the duration of the segment
            segmentDurations = cat(1, segmentDurations, [lastBin - firstBin + 1]);
            
            % Extract EEG data for this segment from this channel
            eegData = EEG.data(chanInd, firstBin:lastBin);
            
            % Extract power
            powerVal = single(multienergyvec(eegData, frequencies, EEG.srate, width));
            
            % Add it to the vector that holds the values for all segments
            % and EEG files
            eegPower = cat(2, eegPower, powerVal);
            
        end % thisSegment
        
        % Add the power values to the cell array
        powerHolder{thisEEG} = eegPower;
        
        % Add the durations to the cell array
        durationHolder{thisEEG} = segmentDurations;
        
    end % thisEEG
    
    % Concatenate the power values from all EEGs
    allPower = [];
    for thisEEG = 1:length(eegList)
        allPower = cat(2, allPower, powerHolder{thisEEG});
    end
    
    % If there's zero power, this will throw a -Inf when you take the log,
    % so remove them
    allPower(allPower == 0) = NaN;
    
    % Take the log of the power values
    logPowerVal = log10(double(allPower));
    
    % Get the mean log(Power) for each frequency
    logPowerMean = nanmean(logPowerVal, 2);
    
    % Fit a chi-square distribution to the power values across
    % frequencies and return the cumulative distribution function
    % (cdf). The first value of the cdf are the frequencies, and
    % the remaining 1000 values are the distribution.
    [cdf, ~] = chi_squarefit(frequencies, logPowerMean);
    cdf = cdf';
    
    % Set the threshold of the distribution
    powerThresh = cdf(:, 10 * amplitudeThresh + 1);
    
    % Loop through each frequency for each segment and return a binary vector that
    % indicates whether the power at that time point exceeded both the 
    % amplitude and duration thresholds
    pepisodeVectors = cell(size(eegList));
    for thisEEG = 1:length(eegList)
        startInd = 1;
        oscillateVector = [];
        
        allPower = powerHolder{thisEEG};
        segmentDurations = durationHolder{thisEEG};
        
        for thisSegment = 1:length(segmentDurations)
            
            % Extract the power data for this segment
            powerData = allPower(:, startInd:startInd + segmentDurations(thisSegment) - 1);
            
            % Update ind
            startInd = startInd + segmentDurations(thisSegment);
            clear tempOscillateVector;
            for thisFreq = 1:length(frequencies)
                tempOscillateVector(thisFreq, :) = single(episodeid(powerData(thisFreq,:), powerThresh(thisFreq), durationThresh * EEG.srate / frequencies(thisFreq), shoulder, 0));
            end % thisFreq
            
            % Add this data to the previous segments
            oscillateVector = cat(2, oscillateVector, tempOscillateVector);
            
        end % thisSegment
        
        % Add this data to the cell array
        oscillateVector = logical(oscillateVector);
        pepisodeVectors{thisEEG} = oscillateVector;
        
    end % thisEEG
    
    % Save the output
    samplerate = EEG.srate;
    save([saveStem chanNames{thisChan} '.mat'], 'pepisodeVectors', 'frequencies','samplerate');
    
end % thisChan



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
    if size(detected,2) < shoulder
        % if our segment is too short, set it to zero
        detected(:) = 0;
    else
        detected(1:shoulder)=0;
        detected((length(detected)-shoulder):length(detected))=0;
    end
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
