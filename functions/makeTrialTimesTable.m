function [trialTimesTable, trialTypeList] = makeTrialTimesTable(EEG, timesNT, timesFT, eventLabelNT, eventLabelFT)
% function trialTimesTable = makeTrialTimesTable(EEG, timesNT, timesFT, eventLabelNT, eventLabelFT)
%
% Make a table that contains for each trial the list of onsetTimes and
% offsetTimes of events of interest. These timings are used to sub-epoch
% each EEG trial. For example, we may want to look at the time interval
% subjects were in the teleporter as well as the second before and the
% second after. Each of these three intervals would be a subepoch defined
% by timesNT and timesFT.
%
% INPUTS:
%   EEG: EEG structure from EEGLAB that contains events
%   timesNT: two-column vector that indictes the onsetTime (col 1) and
%       offsetTime (col 2) of each NT event
%   timesFT: two-column vector that indicates the onsetTime (col 1) and
%       offsetTime (col 2) of each FT event
%   eventLabelNT: the label that corresponds to "NT" events in the EEG
%       structure; in this experiment the time label is always the second
%       digit of the label (e.g., set eventLabelNT = '1' will identify both
%       NSNT trials ('11') and FSNT trials ('21')
%   eventLabelFT: the label that corresponds to "FT" events in the EEG
%       structure; in this experiment the time label is always the second
%       digit of the label (e.g., set eventLabelFT = '2' will identify both
%       NSFT trials ('12') and FSFT trials ('22')
%
% OUTPUTS:
%   trialTimesTable: table with three columns: trial time type, onsetTimes
%       (vector) and offsetTimes (vector)
%   trialTypeList: list of the full trial type (space + time) for each
%       event (e.g., '11')

% Keep only the 2nd value of each trial type since this indicates the time
% type (1 = NT, 2 = FT)
trialTypeList = {EEG.event.type}';
boundaryInds  = find(strcmpi('boundary', trialTypeList));
trialTypeList(boundaryInds) = [];
trialTypeList = cellstr(cellfun(@(s) s(2), trialTypeList));

if (isempty(trialTypeList))
    error('No trial events found.')
end

% Make a table of trial type list
trialTypeTable = table(trialTypeList);

% Make tables of onsets and offsets
onsetTimeTable  = table({eventLabelNT; eventLabelFT}, [timesNT(:, 1)'; timesFT(:, 1)'], 'VariableNames', {'trialTypeList', 'onsetTimes'});
offsetTimeTable = table({eventLabelNT; eventLabelFT}, [timesNT(:, 2)'; timesFT(:, 2)'], 'VariableNames', {'trialTypeList', 'offsetTimes'});

% Join the tables together
trialTimesTable = join(trialTypeTable, onsetTimeTable);
trialTimesTable = join(trialTimesTable, offsetTimeTable);

% Return the original trialTypeList
trialTypeList = {EEG.event.type}';

end