% This script will run the function "Teleporter_Epoched_Pepisode" for each
% session for each subject.
%
% Lindsay Vass
% 2 June 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc;

%% Load structure with session info
sessionInfoPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo.mat';
load(sessionInfoPath);

%% Set up analysis parameters

analysisDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Timing_Onset_Offset/';

% time periods of interest in ms relative to teleporter entry
timePointNames = {'WholeTrial'};

timesNT = [-3000 4830];
timesFT = [-3000 5830];

% frequencies to use
frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% Check that all frequencies are valid for pepisode
intervalsNT = timesNT(:,2) - timesNT(:, 1) + 1;
intervalsFT = timesFT(:,2) - timesFT(:, 1) + 1;
minInterval = min([intervalsNT; intervalsFT]);

durationThresh = 3; % duration for pepisode in cycles
minFrequency = durationThresh / (minInterval / 1000);

excludedFreqs = frequencies(frequencies < minFrequency);
frequencies   = frequencies(frequencies >= minFrequency);

if ~isempty(excludedFreqs)
    fprintf(['\n\n\n WARNING: \n\n']);
    fprintf(['Pepisode cannot be calculated for some frequencies because there \n' ...
        'must be ' num2str(durationThresh) ' complete cycles to estimate pepisode.\n' ...
        'Excluding the following frequencies:\n']);
    fprintf([num2str(excludedFreqs) '\n\n']);
end
%% Run the analysis for each session

for thisSubject = 1:size(sessionInfo, 2)
    
    fprintf('\n\n-------------------------------------------------------\n\n');
    fprintf([sessionInfo(thisSubject).subjectID '\n']);
    
    for thisSession = 1:size(sessionInfo(thisSubject).teleporter, 2)
        
        % Extract session-specific data from structure
        subjectID  = sessionInfo(thisSubject).subjectID;
        teleporter = sessionInfo(thisSubject).teleporter{thisSession};
        chanList   = sessionInfo(thisSubject).chanList;
        
        fprintf([teleporter '\n\n']);
        
        % Set subject directory
        subjectDir = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/'];
        
        % Specify file naming conventions for data. There are
        % separate EEG files for each depth electrode, so we will specify every
        % part of the path except the depth electrode name. Then, we will later
        % combine these as [prefix depthName suffix]. We will use a cell array to
        % allow for multiple prefixes or suffixes (e.g., different prefixes for
        % EDF1 and EDF2)
        numEDFs = sessionInfo(thisSubject).numEDFs(thisSession);
        if numEDFs == 1
            cleanedUnepochedPrefix = {[subjectDir 'PreProcessing Intermediates/' subjectID '_' teleporter '_unepoched_']};
            cleanedUnepochedSuffix = {'_noSpikes_noWaves.set'};
        else
            cleanedUnepochedPrefix = cell(numEDFs, 1);
            for thisEDF = 1:numEDFs
                
                cleanedUnepochedPrefix(thisEDF) = {[subjectDir 'PreProcessing Intermediates/' subjectID '_' teleporter '_EDF' num2str(thisEDF) '_unepoched_']};
                cleanedUnepochedSuffix(thisEDF) = {'_noSpikes_noWaves.set'};
                
            end % thisEDF
        end
        
        
        cleanedEpochedPrefix = {[subjectDir 'Epoched Data/' subjectID '_' teleporter '_epoched_']};
        cleanedEpochedSuffix = {'_noSpikes_noWaves.set'};
        
        
        % Specify path to save the cell array of power values to
        saveFile = [subjectID '_' teleporter '_pepisode_epoched'];
        
        % Run the analysis
        Teleporter_Pepisode_Timing_Onset_Offset(subjectID, ...
            subjectDir, ...
            teleporter, ...
            chanList, ...
            cleanedUnepochedPrefix, ...
            cleanedUnepochedSuffix, ...
            cleanedEpochedPrefix, ...
            cleanedEpochedSuffix, ...
            analysisDir, ...
            saveFile, ...
            timePointNames, ...
            timesNT, ...
            timesFT, ...
            frequencies);
        
    end % thisSession
    
end % thisSubject

