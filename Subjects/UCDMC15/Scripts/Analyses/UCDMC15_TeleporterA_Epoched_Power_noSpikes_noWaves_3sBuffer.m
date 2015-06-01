% This script will extract spectral power for a given set of time points
% and frequencies in seven epochs:
%   - Pre-teleportation (-3000 : -2001 ms) (-2000 : -1001 ms) (-1000 : 0 ms)
%   - Teleportation (0 : 1830 ms for NT or 0 : 2830 ms for FT)
%   - Post-teleportation
%       - NT: (1831 : 2830 ms) (2831 : 3830 ms) (3831 : 4830 ms)
%       - FT: (2831 : 3830 ms) (3831 : 4830 ms) (4831 : 5830 ms)
%
% Lindsay Vass 29 May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% set parameters for analysis

% Subject info
subjectID  = 'UCDMC15';
subjectDir = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/'];
teleporter = 'TeleporterA';

% Specify file naming conventions for data. There are
% separate EEG files for each depth electrode, so we will specify every
% part of the path except the depth electrode name. Then, we will later
% combine these as [prefix depthName suffix]. We will use a cell array to
% allow for multiple prefixes or suffixes (e.g., different prefixes for
% EDF1 and EDF2)
cleanedEpochedPrefix = {[subjectDir 'Epoched Data/' subjectID '_' teleporter '_epoched_']};
cleanedEpochedSuffix = {'_noSpikes_noWaves.set'};

% Specify path to save the power calculations to
saveStem = [subjectDir 'Mat Files/Power/' subjectID '_' teleporter '_epoched_power_'];

% Specify path to save the cell array of power values to
saveFile = [subjectDir 'Mat Files/Power/Summary/' subjectID '_' teleporter '_epoched_power_summary_noSpikes_noWaves_3sBuffer'];

% channel names to use
chanList  = {'RAD3' 'RAD4' 'RAD5' 'RAD6' 'RHD1' 'RHD2' 'RHD3' 'RHD4' 'LAD3' 'LAD4' 'LAD5' 'LHD1' 'LHD2' 'LHD3'};

% time periods of interest in ms relative to teleporter entry
timePointNames = {'Pre3' 'Pre2' 'Pre1' 'Tele' 'Post1' 'Post2' 'Post3'};

timesNT = [-3000 -2001; -2000 -1001; -1000 0; 1 1830; 1831 2830; 2831 3830; 3831 4830];
timesFT = [-3000 -2001; -2000 -1001; -1000 0; 1 2830; 2831 3830; 3831 4830; 4831 5830];

% frequencies to use
frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% set paths and filenames
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/PepisodeCode/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/arne_code/'));

% Make the output directory if it doesn't already exist
if ~exist([subjectDir 'Mat Files/Power/'], 'dir')
    system(['mkdir ' subjectDir 'Mat\ Files/Power/']);
end

if ~exist([subjectDir 'Mat Files/Power/Summary/'], 'dir')
    system(['mkdir ' subjectDir 'Mat\ Files/Power/Summary/']);
end

% Get depth names from chanList
depthNames = unique(cellfun(@(s) s(1:3), chanList, 'UniformOutput', false));

%% Calculate power for each epoch of each trial

for thisDepth = 1:length(depthNames)
    
    % Load EEG data
    eegPath = [cleanedEpochedPrefix{1} depthNames{thisDepth} cleanedEpochedSuffix{1}];
    EEG = pop_loadset(eegPath);
    
    %% Create tables of epoch onsets and offsets
    
    % Keep only the 2nd value of each trial type since this indicates the time
    % type (1 = NT, 2 = FT)
    trialTypeList = {EEG.event.type}';
    trialTypeList = cellstr(cellfun(@(s) s(2), trialTypeList));
    
    % Make a table of trial type list
    trialTypeTable = table(trialTypeList);
    
    % Make tables of onsets and offsets
    onsetTimeTable  = table({'1'; '2'}, [timesNT(:, 1)'; timesFT(:, 1)'], 'VariableNames', {'trialTypeList', 'onsetTimes'});
    offsetTimeTable = table({'1'; '2'}, [timesNT(:, 2)'; timesFT(:, 2)'], 'VariableNames', {'trialTypeList', 'offsetTimes'});
    
    % Join the tables together
    trialTimesTable = join(trialTypeTable, onsetTimeTable);
    trialTimesTable = join(trialTimesTable, offsetTimeTable);
    
    
    % Find the channels on this depth electrode
    chanDepth = char(chanList);
    chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
    chanDepth = cellstr(chanDepth);
    chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
    chanNames = chanList(chanInd);
    
    
    %% Calculate power
    calcEpochedPowerLKV(EEG, chanNames, {EEG.event.type}, trialTimesTable{:, 2}, trialTimesTable{:, 3}, saveStem, frequencies, 6, 3);
    
    
end % thisDepth


%% Extract power values for our epochs of interest

% Initialize the cell array to hold all of our pepisode values
powerSummary      = cell(1,10);
powerSummary(1,:) = {'SubjectID','Teleporter','Electrode','TrialNumber','TrialSpaceType','TrialTimeType','TrialType','TimePoint','Frequency','Power'};
thisRow = 2;

for thisDepth = 1:length(depthNames)
    
    fprintf(['\n\nWorking on electrode #' num2str(thisDepth) ' of ' num2str(length(depthNames)) '\n']);
    
    % Find the channels on this depth electrode
    chanDepth = char(chanList);
    chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
    chanDepth = cellstr(chanDepth);
    chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
    chanNames = chanList(chanInd);
    
    % Loop through channels on this depth electrode
    for thisChan = 1:length(chanNames)
        
        % Path to the power vector we calculated in the previous step
        powerVectorFile = [saveStem chanNames{thisChan} '.mat'];
        
        % Load the power vector file
        load(powerVectorFile);
        
        % Loop through trials
        for thisTrial = 1:size(trialTypes, 1)
            
            % Extract the trial type for this trial
            thisLabel = trialTypes{thisTrial};
            
            % Determine whether it's NS or FS
            if strcmpi('1', thisLabel(1)) == 1
                thisSpaceType = 'NS';
            elseif strcmpi('2', thisLabel(1)) == 1
                thisSpaceType = 'FS';
            else
                error('Unknown trial type')
            end
            
            % Determine whether it's NT or FT
            if strcmpi('1', thisLabel(2)) == 1
                thisTimeType = 'NT';
            elseif strcmpi('2', thisLabel(2)) == 1
                thisTimeType = 'FT';
            else
                error('Unknown trial type')
            end
            
            % Combine them together to make the spatiotemporal type
            thisType = [thisSpaceType thisTimeType];
            
            % Loop through time points
            for thisTimePoint = 1:length(timePointNames)
                
                powerVec = allPowerData{thisTrial, thisTimePoint};
                
                % Take the mean across time
                meanPower = nanmean(powerVec, 2);
                
                % Add the values to the summary cell array
                for thisFreq = 1:length(frequencies)
                    
                    powerSummary(thisRow,:) = {subjectID, teleporter, chanNames{thisChan}, thisTrial, thisSpaceType, thisTimeType, thisType, timePointNames{thisTimePoint}, frequencies(thisFreq), meanPower(thisFreq)};
                    
                    thisRow = thisRow + 1;
                    
                end % thisFreq
                
            end % thisTimePoint
            
        end % thisTrial
        
    end % thisChan
    
    
end % thisDepth

% Write out the summary cell array to file
fprintf('\n\nWriting output to file...\n');
dlmcell([saveFile '.csv'], powerSummary, 'delimiter', ',');
save([saveFile '.mat'], 'powerSummary');
