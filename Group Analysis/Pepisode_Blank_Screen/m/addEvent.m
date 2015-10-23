function EEG = addEvent(EEG, onset, label)
% Add event to the EEG structure
% >> EEG = addEvent(EEG, onset, label)
%
% Inputs:
%   EEG: EEGLAB structure
%   onset: time of the event onset in EEG samples
%   label: numeric value used to label the event
%
% Output:
%   EEG: EEGLAB structure with added event

totalEvents = size(EEG.event, 2);
EEG.event(totalEvents + 1).latency = onset;
EEG.event(totalEvents + 1).type = label;