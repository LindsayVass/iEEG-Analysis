function [boundaries] = findBoundaries(EEG)

% function [boundaries] = findBoundaries(EEG)
% Extract the boundary events from the EEG structure
%
% INPUT: 
%   EEG: EEG structure from EEGLAB
%
% OUTPUT:
%   boundaries: vector of event onsets in EEG samples

% get the names of all the events in this dataset
eventLabels = {EEG.event.type};

% find the boundary events
boundaryInds = strcmpi('boundary', eventLabels);

% Get the onsets in EEG samples of the boundary events
boundaries = floor(cell2mat({EEG.event(boundaryInds).latency}));

end