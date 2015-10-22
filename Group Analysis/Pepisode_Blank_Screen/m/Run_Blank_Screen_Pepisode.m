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
        
        % Loop over depth electrodes
        for thisElec = 1:length(sessionInfo(thisSubject).teleporter(thisSession).depths)
            
            % Make epochs
            
            % Calculate pepisode
            
            % Export data
            
        end % thisElec
        
    end % thisSession
    
end % thisSubject