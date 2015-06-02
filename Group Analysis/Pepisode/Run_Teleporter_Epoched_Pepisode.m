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
% time periods of interest in ms relative to teleporter entry
timePointNames = {'Pre3' 'Pre2' 'Pre1' 'Tele' 'Post1' 'Post2' 'Post3'};

timesNT = [-3000 -2001; -2000 -1001; -1000 0; 1 1830; 1831 2830; 2831 3830; 3831 4830];
timesFT = [-3000 -2001; -2000 -1001; -1000 0; 1 2830; 2831 3830; 3831 4830; 4831 5830];

% frequencies to use
frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% Run the analysis for each session

for thisSubject = 1:size(sessionInfo, 2)
    
    for thisSession = 1:size(sessionInfo(thisSubject).teleporter, 2)
        
        % Extract session-specific data from structure
        subjectID  = sessionInfo(thisSubject).subjectID;
        teleporter = sessionInfo(thisSubject).teleporter{thisSession};
        chanList   = sessionInfo(thisSubject).chanList;
        
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
        saveFile = [subjectDir 'Mat Files/Pepisode/Summary/' subjectID '_' teleporter '_epoched_pepisode_summary_noSpikes_noWaves_3sBuffer'];
        
        % Run the analysis
        Teleporter_Epoched_Pepisode(subjectID, ...
            subjectDir, ...
            teleporter, ...
            chanList, ...
            cleanedUnepochedPrefix, ...
            cleanedUnepochedSuffix, ...
            cleanedEpochedPrefix, ...
            cleanedEpochedSuffix, ...
            saveFile, ...
            timePointNames, ...
            timesNT, ...
            timesFT, ...
            frequencies);
        
    end % thisSession
    
end % thisSubject

