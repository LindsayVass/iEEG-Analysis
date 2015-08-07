function epochOnsets = getEpochOnsetsPTB(timeSyncPath, systemTimeStart)
% epochOnsets = getEpochOnsetsPTB(timeSyncPath, systemTimeStart)
%
% Purpose: Take the timing of epochs onsets in Unix ticks and convert them
%   to indices of the EEG data.
%
% INPUT:
%   timeSyncPath: path to the .mat file containing the parameters for the
%       linear regression relating ticks to EEG bins
%   systemTimeStart: vector of epoch onset times in Unix ticks; output from
%       "sampleBehavioralData.m"
%
% OUTPUT:
%   epochOnsets: vector of onset times in EEG indices
%
% Author: Lindsay Vass
% Date: 5 August 2015

load(timeSyncPath);

epochOnsets = nan(size(systemTimeStart));

for thisInterval = 1:length(epochOnsets)
    epochOnsets(thisInterval) = round(systemTimeStart(thisInterval) * time_sync_regression(1) + time_sync_regression(2));
end

end

