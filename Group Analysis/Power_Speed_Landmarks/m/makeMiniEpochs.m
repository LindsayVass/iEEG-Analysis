function [newEEG, goodEpochs] = makeMiniEpochs(epochOnsets, epochLabel, epochInterval, unepochedEEGPath)
% [newEEG, goodEpochs] = makeMiniEpochs(epochOnsets, epochLabel, epochInterval, unepochedEEGPath)
% 
% Purpose: Insert epoch events into an EEGLAB structure that has been
%          marked for artifacts and then epoch the EEG data.
%
% INPUT:
%   epochOnsets: a vector of epoch start times in EEG bins; output of
%       either "getEpochOnsets" or "getEpochOnsetsPTB"
%   epochLabel: a numeric or string used to label the epochs inserted into
%       the data
%   epochInterval: two-value vector indicating the start and stop times of
%       the epoch in seconds, relative to the onset times in epochOnsets;
%       for example, [0 0.2] is a 200 ms epoch that starts at the times
%       indicated in epochOnsets and ends 200 ms later
%   unepochedEEGPath: path to the .set file that has been previously marked
%       during artifact detection
%
% OUTPUT:
%   newEEG: EEGLAB structure containing the epoched data
%   goodEpochs: indices of the epochs in epochOnsets which were
%       subsequently included in newEEG (i.e., did not overlap with an
%       artifact)
%
% Author: Lindsay Vass
% Date: 5 August 2015

EEG = pop_loadset(unepochedEEGPath);
EEG = addEventsToEEG(EEG, epochOnsets, epochLabel);
[EEG, noWaveEpochs]  = removeArtifactsFromEEG(EEG);
[newEEG, goodEpochs] = pop_epoch(EEG, {epochLabel}, epochInterval);
goodEpochs = noWaveEpochs(goodEpochs);

end

function EEG = addEventsToEEG(EEG, epochOnsets, label)

% add epochs to the EEG event struct
thisEpoch = 1;
numSpikes = size(EEG.event, 2);
for n = numSpikes + 1:numSpikes + length(epochOnsets)
    EEG.event(n).latency = epochOnsets(thisEpoch);
    EEG.event(n).type = label;
    thisEpoch = thisEpoch + 1;
end
end

function [noWaveEEG, noWaveEpochs] = removeArtifactsFromEEG(EEG)

% remove spikes and waves
artifactLabels = {'spike', 'complex', 'other', 'sharpWave'};
indexes = marks_label2index(EEG.marks.time_info, artifactLabels, 'indexes', 'exact', 'on');
origEventLatency = cell2mat({EEG.event.latency});

if isempty(indexes)
    noWaveEEG = EEG;
    warning('No sharp waves found.')
else
    [noWaveEEG, ~] = pop_marks_select_data(EEG, 'time marks', indexes);
end

% if the spikes or waves overlap with a real event, it will be removed,
% so check for that now
[~, missingWaveInd] = intersect(origEventLatency, indexes);
noWaveEpochs = [1:1:length(origEventLatency)]';
noWaveEpochs(missingWaveInd) = [];

end