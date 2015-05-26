function powerVec = getPowerLKV(powerFile,EDFnumber,latencyBin,durationMS,offsetMS,frequencies)
% getPowerLKV - Return the zscored power data from a file for a set
% of desired events. 
%
% FUNCTION powerVec=getPowerLKV(powerFile, EDFnumber, latencyBin, durationMS, offsetMS,frequencies)
%
% INPUT ARGs:
%   powerFile - the file calculated by calcPowerLKV.m with the
%       data for this set of events
%   EDFnumber - default = 1; which EDF we're evaluating
%   latencyBin - time in EEG samples of the event of interest
%   durationMS - length of signal to extract in milliseconds
%   offsetMS - offset in milliseconds for the start time, relative to
%       the latency
%   frequencies - vector of frequencies
%
%
% OUTPUT ARGs:
%   powerVec - (Events,Freqs,Time) - vector of zscored power at that point 
%       in time
%


% Get the data and determine durations in EEG samples
load(powerFile);
duration = round((durationMS)*samplerate/1000);
offset = round((offsetMS)*samplerate/1000);

% Select the appropriate data based on the EDFnumber
powerVector = powerVectors{EDFnumber};

% Allocate the output vector
powerVec = nan(length(latencyBin),length(frequencies),duration);

for thisEvent = 1:length(latencyBin)
    
    % find the right point in the array of pepisode values
    startTime = offset + latencyBin(thisEvent);
    
    if startTime > 0
        if startTime + duration - 1 > length(powerVector)
            warning('Dataset does not include requested time points. Returning NaNs.');
            continue
        else
            powerVec(thisEvent,:,:)= powerVector(:,startTime:(startTime+duration-1));
        end
    else
        error('Invalid index. %n + %n is not > 0.',offset,latencyBin(thisEvent));
    end
end
    
