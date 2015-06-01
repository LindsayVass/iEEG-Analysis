% This script will extract spectral power for a given set of time points 
% and frequencies in seven epochs:
%   - Pre-teleportation (-3000 : -2001 ms) (-2000 : -1001 ms) (-1000 : 0 ms)
%   - Teleportation (0 : 1830 ms for NT or 0 : 2830 ms for FT)
%   - Post-teleportation 
%       - NT: (1831 : 2830 ms) (2831 : 3830 ms) (3831 : 4830 ms)
%       - FT: (2831 : 3830 ms) (3831 : 4830 ms) (4831 : 5830 ms)
%
% Lindsay Vass 26 May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% set parameters for analysis

% Subject info
subjectID  = 'UCDMC14';
subjectDir = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/'];
teleporter = 'TeleporterA';

% Specify file naming conventions for data. There are
% separate EEG files for each depth electrode, so we will specify every
% part of the path except the depth electrode name. Then, we will later
% combine these as [prefix depthName suffix]. We will use a cell array to
% allow for multiple prefixes or suffixes (e.g., different prefixes for
% EDF1 and EDF2)
cleanedUnepochedPrefix = {[subjectDir 'PreProcessing Intermediates/' subjectID '_' teleporter '_unepoched_']};
cleanedUnepochedSuffix = {'_noSpikes_noWaves.set'};

% Specify path to save the power calculations to
saveStem = [subjectDir 'Mat Files/Power/' subjectID '_' teleporter '_power_'];

% Specify path to save the cell array of power values to
saveFile = [subjectDir 'Mat Files/Power/Summary/' subjectID '_' teleporter '_power_summary_noSpikes_noWaves_3sBuffer'];

% channel names to use
chanList  = {'LAD1' 'LHD1' 'LHD2' 'RAD1' 'RHD1' 'RHD2'};

% if the power across the recording has already been
% calculated, set this to 1
skipCompute = 0;

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


%% Calculate power distributions

% In this first step, we will use the cleaned, unepoched EEG files and
% calculate the z-score of the log(power) of the EEG data at each
% frequency.

if ~skipCompute
    
    for thisDepth = 1:length(depthNames)
        
        % Find the channels on this depth electrode
        chanDepth = char(chanList);
        chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
        chanDepth = cellstr(chanDepth);
        chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
        chanNames = chanList(chanInd);
        
        % Set up eeg list for this channel
        numEDFs = size(cleanedUnepochedPrefix,2);
        eegList = {};
        for thisEDF = 1:numEDFs
            pathName = [cleanedUnepochedPrefix{thisEDF} depthNames{thisDepth} cleanedUnepochedSuffix{1}];
            eegList{thisEDF} = pathName;
            
        end % thisEDF
        
        % Calculate pepisode
        calcPowerLKV(eegList, chanNames, saveStem, frequencies, 6);
        
        
    end % thisDepth
    
end

%% Extract power values for our epochs of interest

% Initialize the cell array to hold all of our pepisode values
powerSummary      = cell(1,11);
powerSummary(1,:) = {'SubjectID','Teleporter','EDF','Electrode','TrialNumber','TrialSpaceType','TrialTimeType','TrialType','TimePoint','Frequency','Power'};
thisRow = 2;

for thisDepth = 1:length(depthNames)
    
    fprintf(['\n\nWorking on electrode #' num2str(thisDepth) ' of ' num2str(length(depthNames)) '\n']);
    
    % Find the channels on this depth electrode
    chanDepth = char(chanList);
    chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
    chanDepth = cellstr(chanDepth);
    chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
    chanNames = chanList(chanInd);
    
    numEDFs = size(cleanedUnepochedPrefix,2);
    
    % Loop through EEG files
    for thisEDF = 1:numEDFs
        
        % Load cleaned unepoched data
        EEG = pop_loadset([cleanedUnepochedPrefix{thisEDF} depthNames{thisDepth} cleanedUnepochedSuffix{1}]);
        
        % Extract the trial events from the event list (i.e., remove boundary
        % events)
        boundaryInds = strcmpi('boundary', {EEG.event.type});
        trialList    = EEG.event(boundaryInds == 0); % select events that are NOT boundaries
        
        % Loop through channels on this depth electrode
        for thisChan = 1:length(chanNames)
            
            % Path to the power vector we calculated in the previous step
            powerVectorFile = [saveStem chanNames{thisChan} '.mat'];
            
            % Load the power vector file
            load(powerVectorFile);
            
            % Loop through trials
            for thisTrial = 1:size(trialList, 2)
                
                % Extract the trial type for this trial
                thisLabel = trialList(thisTrial).type;
                
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
                    
                    if strcmpi('NT', thisTimeType) == 1
                        
                        powerVec = getPowerLKV(powerVectors, samplerate, thisEDF, trialList(thisTrial).latency, timesNT(thisTimePoint, 2) - timesNT(thisTimePoint, 1) , timesNT(thisTimePoint, 1), frequencies);
                        
                    else
                        
                        powerVec = getPowerLKV(powerVectors, samplerate, thisEDF, trialList(thisTrial).latency, timesFT(thisTimePoint, 2) - timesFT(thisTimePoint, 1) , timesFT(thisTimePoint, 1), frequencies);
                        
                    end
                    
                    % Take the mean across time
                    meanPower = nanmean(powerVec,3);
                    
                    % Add the values to the summary cell array
                    for thisFreq = 1:length(frequencies)
                        
                        powerSummary(thisRow,:) = {subjectID, teleporter, thisEDF, chanNames{thisChan}, thisTrial, thisSpaceType, thisTimeType, thisType, timePointNames{thisTimePoint}, frequencies(thisFreq), meanPower(thisFreq)};
                        
                        thisRow = thisRow + 1;
                        
                    end % thisFreq
                    
                end % thisTimePoint
                 
            end % thisTrial
            
        end % thisChan
        
    end % thisEDF
    
end % thisDepth

% Write out the summary cell array to file
fprintf('\n\nWriting output to file...\n');
dlmcell([saveFile '.csv'], powerSummary, 'delimiter', ',');
save([saveFile '.mat'], 'powerSummary');
