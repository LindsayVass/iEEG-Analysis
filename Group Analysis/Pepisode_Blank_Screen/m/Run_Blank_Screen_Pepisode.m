% The purpose of this script is to measure pepisode in the period
% immediately after free exploration and navigation when participants are
% looking at a black screen containing some text. This script will first
% identify the onset times from the behavioral text files, convert those
% from ticks to EEG samples, epoch the data, calculate pepisode, and then
% export the data.
%
% Lindsay Vass
% 22 October 2015


% Structure containing subject/session info
experimentPath  = '/Users/Lindsay/Documents/MATLAB/iEEG/';
sessionInfoPath = [experimentPath 'Group Analysis/Subject Info/SessionInfo2.mat'];
load(sessionInfoPath);

% Loop over subjects
for thisSubject = 1:length(sessionInfo.subjectID)
    
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
            
            
            % Make epochs
            
            % Calculate pepisode
            
            % Export data
            
        end % thisElec
        
    end % thisSession
    
end % thisSubject