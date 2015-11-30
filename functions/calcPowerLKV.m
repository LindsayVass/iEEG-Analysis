function allPower = calcPowerLKV(eegList, chanName, frequencies, waveletWidth)

% Return power at each time point at a given set of
% frequencies using Morlet wavelets.
%
% >> allPower = calcPowerLKV(eegList, chanName, frequencies, waveletWidth, saveFile)
%
% INPUTS:
%   eegList: cell array of paths to the EEGLAB dataset(s) on which to
%       calculate power
%   chanName: string of electrode name to use for analysis
%   frequencies: vector of frequencies in Hz at which to calculate power
%       distributions
%   waveletWidth (optional): number of Morlet wavelet cycles to use for
%       power estimation; default = 6
%
% OUTPUT:
%   allPower: frequencies x time points matrix of power values
%
% Lindsay Vass

% Specify waveletWidth if not indicated by user
if nargin < 4
    warning('Argument waveletWidth not specified. Using default value 6.');
    waveletWidth = 6;
end

if nargin < 3
    error('Arguments eegList, chanNames, and frequencies are required.')
end

% Make sure the EEG files exist
for thisFile = 1:length(eegList)
    if ~exist(eegList{thisFile}, 'file')
        error(['EEG file does not exist: ' eegList{thisFile}])
    end
end

% Initialize output struct
powerDistribution = struct;

% Initialize struct to hold all of the power values from each segment of
% each EEG
powerHolder = struct;

% Loop through EEG datasets
for thisEEG = 1:length(eegList)
    
    % Load the EEG data
    EEG = pop_loadset(eegList{thisEEG});
    
    % Identify the boundaries within the EEG file
    boundaries = findBoundaries(EEG);
    
    % find the index of this channel
    chanInd = find(strcmpi(chanName, {EEG.chanlocs.labels}));
    
    % Calculate power
    [eegPower] = extractPowerWithBoundaries(EEG, chanInd, boundaries, frequencies, waveletWidth);
    
    % Add to summary array
    powerHolder(thisEEG).EEG = eegPower;
    
end % thisEEG

% Build the power distribution by concatenating all power values
allPower = [];

for thisEEG = 1:size(powerHolder, 2)
    
    for thisSegment = 1:size(powerHolder(thisEEG).EEG, 2)
        
        allPower = cat(2, allPower, powerHolder(thisEEG).EEG(thisSegment).segment.power);
        
    end % thisSegment
    
end % thisEEG

% If there's zero power, this will cause an error when we take the log, so
% remove them
allPower(allPower == 0) = NaN;

