function calcPowerLKV(eegList, chanNames, saveStem, frequencies, width)
% calcPowerLKV: Calculate power on the clean segments of the EEG.
% This function will loop through the channels in "chanNames" and calculate
% the z-score of the power of the EEG data at each frequency. 
% The function will first loop through all EEG files and calculate the 
% power at each time point. It will then log-transform the power and
% normalize it by z-scoring within each electrode and frequency with
% respect to the mean and standard deviation. It will save these values to 
% an output file where they can be accessed later using getPowerLKV.
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
%       '/path/to/power/UCDMC15_TeleporterA_power_'
%
%   width: number of cycles to use for wavelet convolution; default = 6
%

% set wavelet width if not specified
if nargin == 4
    width = 6;
end

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
    
    % Get the standard deviation for each frequency
    logPowerSD = nanstd(logPowerVal, 0, 2);
    
    % Calculate the zscore of each logPower value
    logPowerZscore = (logPowerVal - repmat(logPowerMean, 1, size(logPowerVal,2)))./repmat(logPowerSD, 1, size(logPowerVal,2));
    
    
    % Loop through each frequency for each segment and return a vector with
    % the zscore of the log(power) at each time point 
    powerVectors = cell(size(eegList));
    for thisEEG = 1:length(eegList)
        startInd = 1;
        powerVector = [];
        
        allPower = powerHolder{thisEEG};
        segmentDurations = durationHolder{thisEEG};
        
        for thisSegment = 1:length(segmentDurations)
            
            % Extract the power zscores for this segment
            powerData = logPowerZscore(:, startInd:startInd + segmentDurations(thisSegment) - 1);
            
            % Update ind
            startInd = startInd + segmentDurations(thisSegment);
            
            % Add this data to the previous segments
            powerVector = cat(2, powerVector, powerData);
            
        end % thisSegment
        
        % Add this data to the cell array
        powerVectors{thisEEG} = powerVector;
        
    end % thisEEG
    
    % Save the output
    samplerate = EEG.srate;
    save([saveStem chanNames{thisChan} '.mat'], 'powerVectors', 'frequencies', 'samplerate');
    
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
