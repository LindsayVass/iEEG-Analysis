function binaryVec = getPepisodeLKV(unionFile,EDFnumber,latencyBin,durationMS,offsetMS,frequencies)
% getPepisodeLKV - Return the Pepisode union data from a file for a set
% of desired events. This function gives the complete vector of zeros
% and ones, you still have to do the averaging to come to a value for
% pepisode (e.g. pepisode=mean(unionvec,3)=0.21)
%
% FUNCTION binaryVec=getPepisodeLKV(unionFile,EDFnumber,latencyBin,durationMS,offsetMS,frequencies)
%
% INPUT ARGs:
%   unionFile - the file calculated by calcPepisodeLKV.m with the
%       data for this set of events
%   EDFnumber - default = 1; which EDF we're evaluating
%   latencyBin - time in EEG samples of the event of interest
%   durationMS - length of signal to extract in milliseconds
%   offsetMS - offset in milliseconds for the start time, relative to
%       the latency
%
%
% OUTPUT ARGs:
%   binaryVec - (Events,Freqs,Time) - vector of 0 and 1 for whether there 
%       is a significant oscillation at that point in time
%


% Get the data and determine durations in EEG samples
load(unionFile);
duration = round((durationMS)*samplerate/1000);
offset = round((offsetMS)*samplerate/1000);

% Select the appropriate data based on the EDFnumber
unionvector = pepisodeVectors{EDFnumber};

% Allocate the output vector
binaryVec=  nan(length(latencyBin),length(frequencies),duration);

for thisEvent = 1:length(latencyBin)
    
    % find the right point in the array of pepisode values
    startTime = offset + latencyBin(thisEvent);
    
    if startTime > 0
        if startTime + duration - 1 > length(unionvector)
            warning('Dataset does not include requested time points. Returning NaNs.');
            continue
        else
            binaryVec(thisEvent,:,:)= unionvector(:,startTime:(startTime+duration-1));
        end
    else
        error('Invalid index. %n + %n is not > 0.',offset,latencyBin(thisEvent));
    end
end
    
