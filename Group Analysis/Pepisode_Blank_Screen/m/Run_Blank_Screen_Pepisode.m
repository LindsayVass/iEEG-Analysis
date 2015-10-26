% The purpose of this script is to measure pepisode in the period
% immediately after free exploration and navigation when participants are
% looking at a black screen containing some text. This script will first
% identify the onset times from the behavioral text files, convert those
% from ticks to EEG samples, epoch the data, calculate pepisode, and then
% export the data.
%
% Lindsay Vass
% 22 October 2015

clear all; close all; clc;

% Structure containing subject/session info
experimentPath  = '/Users/Lindsay/Documents/MATLAB/iEEG/';
sessionInfoPath = [experimentPath 'Group Analysis/Subject Info/SessionInfo2.mat'];
load(sessionInfoPath);

% Epoch intervals in ms
shortEpoch = [0 1.83];
longEpoch  = [0 2.83];

% Which artifacts to exclude
excludedArtifacts = {'spike','complex','other','sharpWave'};

% Prepare output array
output = {'Subject', 'Teleporter', 'Electrode', 'Condition', 'EpochLength', 'Frequency', 'Pepisode'};

% Path to save csv output
outputPath = [experimentPath 'Group Analysis/Pepisode_Blank_Screen/csv/output_' date '.csv'];

% Loop over subjects
for thisSubject = 1:length(sessionInfo)
    
    subjectID = sessionInfo(thisSubject).subjectID{1};
    
    % Loop over sessions
    for thisSession = 1:length(sessionInfo(thisSubject).teleporter)
        
        sessionID = sessionInfo(thisSubject).teleporter(thisSession).name{1};
        
        % Set paths to behavioral text files
        if strcmpi(subjectID, 'UCDMC15') == 0
            
            % handle different naming for each session
            if thisSession == 1
                sessionTxt = '';
            else
                sessionTxt = ' 2';
            end
            
            freeExplorePath = [experimentPath 'Subjects/' subjectID '/Behavioral Data/' sessionID '/s' num2str(thisSubject) '_freeexplore_patientTeleporterData' sessionTxt '.txt'];
            navigationPath  = [experimentPath 'Subjects/' subjectID '/Behavioral Data/' sessionID '/s' num2str(thisSubject) '_patientTeleporterData' sessionTxt '.txt'];
            
        else % handle different file name for UCDMC15
            freeExplorePath = [experimentPath 'Subjects/' subjectID '/Behavioral Data/' sessionID '/s' num2str(thisSubject) '_FreeExplore_' sessionID '.txt'];
            navigationPath  = [experimentPath 'Subjects/' subjectID '/Behavioral Data/' sessionID '/s' num2str(thisSubject) '_FindStore_' sessionID '_FIXED.txt'];
        end
        
        % Parse text files to get onset times for end of free explore and
        % end of navigation
        if strcmpi(subjectID, 'UCDMC15') == 0
            
            freeExploreTick = getFreeExploreTick(freeExplorePath, 0);
            navigationTick  = getNavigationTick(navigationPath, 0);
            
        else % handle different txt file format for UCDMC15
            
            freeExploreTick = getFreeExploreTick(freeExplorePath, 1);
            navigationTick  = getNavigationTick(navigationPath, 1);
            
        end
        
        % Conver time in ticks to time in EEG samples
        if strcmpi(subjectID, 'UCDMC15') == 0
            
            % load time sync mat file
            load([experimentPath 'Subjects/' subjectID '/Mat Files/' subjectID '_' sessionID  '_time_sync.mat'])
            
            % convert from ticks to samples
            freeExploreOnset = round(freeExploreTick * time_sync_regression(1) + time_sync_regression(2));
            navigationOnset  = round(navigationTick * time_sync_regression(1) + time_sync_regression(2));
            
            
        else % use unity pulses for UCDMC15
            
            % load pulses mat files (one for each EDF)
            load([experimentPath 'Subjects/' subjectID '/Mat Files/' subjectID '_' sessionID '_EDF1_pulse_timing.mat'])
            indEEG1 = indEEG;
            unityTicks1 = unityTicks;
            
            load([experimentPath 'Subjects/' subjectID '/Mat Files/' subjectID '_' sessionID '_EDF2_pulse_timing.mat'])
            indEEG2 = indEEG;
            unityTicks2 = unityTicks;
            
            clear indEEG unityTicks;
            
            % convert from ticks to samples
            [freeExploreOnset, freeExploreEDF] = tick2sampleUnity(freeExploreTick, unityTicks1, unityTicks2, indEEG1, indEEG2);
            [navigationOnset, navigationEDF]  = tick2sampleUnity(navigationTick, unityTicks1, unityTicks2, indEEG1, indEEG2);
            
        end
        
        % Loop over depth electrodes
        for thisElec = 1:length(sessionInfo(thisSubject).teleporter(thisSession).depths)
            
            clear freeExploreEEG navigationEEG shortFreeExploreEEG longFreeExploreEEG shortNavigationEEG longNavigationEEG
            
            depthID  = sessionInfo(thisSubject).teleporter(thisSession).depths(thisElec).name{1};
            chanList = sessionInfo(thisSubject).teleporter(thisSession).depths(thisElec).chanList;
            
            % get path of EEG file
            if strcmpi(subjectID, 'UCDMC15') == 0
                freeExploreEegPath = [experimentPath 'Subjects/' subjectID '/PreProcessing Intermediates/' subjectID '_' sessionID '_unepoched_' depthID '_marked.set'];
                navigationEegPath  = [experimentPath 'Subjects/' subjectID '/PreProcessing Intermediates/' subjectID '_' sessionID '_unepoched_' depthID '_marked.set'];
            else % handle 2 EDFs
                freeExploreEegPath = [experimentPath 'Subjects/' subjectID '/PreProcessing Intermediates/' subjectID '_' sessionID '_EDF' num2str(freeExploreEDF) '_unepoched_' depthID '_marked.set'];
                navigationEegPath = [experimentPath 'Subjects/' subjectID '/PreProcessing Intermediates/' subjectID '_' sessionID '_EDF' num2str(navigationEDF) '_unepoched_' depthID '_marked.set'];
            end
            
            % load EEG
            freeExploreEEG = pop_loadset(freeExploreEegPath);
            navigationEEG  = pop_loadset(navigationEegPath);
            
            % Insert events into EEG
            freeExploreEEG = addEvent(freeExploreEEG, freeExploreOnset, 1);
            freeExploreEEG = addEvent(freeExploreEEG, freeExploreOnset, 2);
            
            navigationEEG  = addEvent(navigationEEG, navigationOnset, 1);
            navigationEEG  = addEvent(navigationEEG, navigationOnset, 2);
            
            % Exclude artifacts
            freeExploreIndexes = marks_label2index(freeExploreEEG.marks.time_info, excludedArtifacts, 'indexes', 'exact', 'on');
            
            if isempty(freeExploreIndexes)
                warning('No artifacts found.')
            else
                [freeExploreEEG,~] = pop_marks_select_data(freeExploreEEG, 'time marks', freeExploreIndexes);
            end
            
            navigationIndexes = marks_label2index(navigationEEG.marks.time_info, excludedArtifacts, 'indexes', 'exact', 'on');
            
            if isempty(navigationIndexes)
                warning('No artifacts found.')
            else
                [navigationEEG,~] = pop_marks_select_data(navigationEEG, 'time marks', navigationIndexes);
            end
            
            % Make epochs
            try
                shortFreeExploreEEG = pop_epoch(freeExploreEEG, {'1'}, shortEpoch);
            catch
                shortFree = 0;
            end
            
            try
                longFreeExploreEEG  = pop_epoch(freeExploreEEG, {'2'}, longEpoch);
            catch
                longFree = 0;
            end
            
            try
                shortNavigationEEG = pop_epoch(navigationEEG, {'1'}, shortEpoch);
            catch
                shortNav = 0;
            end
            
            try
                longNavigationEEG  = pop_epoch(navigationEEG, {'2'}, longEpoch);
            catch
                longNav = 0;
            end
            
            
            % Loop over channels on this electrode
            chanList = sessionInfo(thisSubject).teleporter(thisSession).depths(thisElec).chanList;
            for thisChan = 1:length(chanList)
                
                chanInd = find(strcmpi(chanList{thisChan}, {freeExploreEEG.chanlocs.labels}));
                
                powerDistPath = [experimentPath 'Subjects/' subjectID '/Mat Files/Pepisode/Power Distributions/' subjectID '_' sessionID '_' chanList{thisChan} '_power_distribution.mat'];
                load(powerDistPath);
                
                % Calculate pepisode
                if ~exist('shortFree', 'var')
                    [~, shortFreeExplorePepisode] = calcEpochedPepisodeLKV(powerDistribution, [powerDistribution.frequency], shortFreeExploreEEG.data(chanInd, :), shortFreeExploreEEG.srate, 95, 3);
                end
                if ~exist('longFree', 'var')
                    [~, longFreeExplorePepisode]  = calcEpochedPepisodeLKV(powerDistribution, [powerDistribution.frequency], longFreeExploreEEG.data(chanInd, :), longFreeExploreEEG.srate, 95, 3);
                end
                if ~exist('shortNav', 'var')
                    [~, shortNavigationPepisode]  = calcEpochedPepisodeLKV(powerDistribution, [powerDistribution.frequency], shortNavigationEEG.data(chanInd, :), shortNavigationEEG.srate, 95, 3);
                end
                if ~exist('longNav', 'var')
                    [~, longNavigationPepisode]   = calcEpochedPepisodeLKV(powerDistribution, [powerDistribution.frequency], longNavigationEEG.data(chanInd, :), longNavigationEEG.srate, 95, 3);
                end
                
                % Add data to output array
                % Loop over frequencies
                freqList = [powerDistribution.frequency];
                for thisFreq = 1:length(freqList)
                    if ~exist('shortFree', 'var')
                        output(end + 1, :) = {subjectID, sessionID, chanList{thisChan}, 'FreeExplore', 'Short', freqList(thisFreq), shortFreeExplorePepisode(thisFreq)};
                    end
                    if ~exist('longFree', 'var')
                        output(end + 1, :) = {subjectID, sessionID, chanList{thisChan}, 'FreeExplore', 'Long', freqList(thisFreq), longFreeExplorePepisode(thisFreq)};
                    end
                    if ~exist('shortNav', 'var')
                        output(end + 1, :) = {subjectID, sessionID, chanList{thisChan}, 'Navigation', 'Short', freqList(thisFreq), shortNavigationPepisode(thisFreq)};
                    end
                    if ~exist('longNav', 'var')
                        output(end + 1, :) = {subjectID, sessionID, chanList{thisChan}, 'Navigation', 'Long', freqList(thisFreq), longNavigationPepisode(thisFreq)};
                    end
                    
                end % thisFreq
                
            end % thisChan
            
        end % thisElec
        
    end % thisSession
    
end % thisSubject

% Save output data
dlmcell(outputPath, output, 'delimiter', ',');