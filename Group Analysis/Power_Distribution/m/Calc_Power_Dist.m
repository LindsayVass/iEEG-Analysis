clear all; close all; clc;

%% Load structure with session info
sessionInfoPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo.mat';
load(sessionInfoPath);

%% Set up analysis parameters
% time periods of interest in ms relative to teleporter entry
% timePointNames = {'Pre1' 'Tele' 'Post1'};
%
timesNT = [-1829 0; 1 1830; 1831 3660];
timesFT = [-2829 0; 1 2830; 2831 5660];

% frequencies to use
frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

% where to save the output
analysisDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power_Distribution/';

waveletWidth = 6;

%% Get power
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
        
        % Get depth names from chanList
        depthNames = unique(cellfun(@(s) s(1:3), chanList, 'UniformOutput', false));
        
        % Loop through depths
        for thisDepth = 1:length(depthNames)
            % make a cell array of the EEG dataset paths
            eegList = cell(length(cleanedUnepochedPrefix), 1);
            for thisEEG = 1:length(cleanedUnepochedPrefix)
                eegPath = [cleanedUnepochedPrefix{thisEEG} depthNames{thisDepth} cleanedUnepochedSuffix{thisEEG}];
                eegList(thisEEG) = {eegPath};
            end
            
            % Find the channels on this depth electrode
            chanNames = findChansOnThisElectrode(chanList, depthNames{thisDepth});
            
            % Loop through channels
            for thisChan = 1:length(chanNames)
                saveFile = [analysisDir 'mat/' sessionInfo(thisSubject).subjectID '_' sessionInfo(thisSubject).teleporter{thisSession} '_' chanNames{thisChan} '.mat'];
                allPower = calcPowerLKV(eegList, chanNames{thisChan}, frequencies, waveletWidth);
                save(saveFile, 'allPower', 'frequencies');
            end
        end
    end
end