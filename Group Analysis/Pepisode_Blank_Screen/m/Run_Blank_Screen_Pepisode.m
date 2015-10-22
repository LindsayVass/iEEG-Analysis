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
sessionInfoPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo2.mat';
load(sessionInfoPath);

% Loop over subjects
for thisSubject = 1:length(sessionInfo.subjectID)
    
    subjectID = sessionInfo(thisSubject).subjectID{1};
    
    % Loop over sessions
    for thisSession = 1:length(sessionInfo(thisSubject).teleporter)
        
        sessionID = sessionInfo(thisSubject).teleporter(thisSession).name{1};
        
        % Set paths to behavioral text files
        
        
        % Loop over depth electrodes
        for thisElec = 1:length(sessionInfo(thisSubject).teleporter(thisSession).depths)
            
            % Get onset times for end of free explore and end of navigation
            
            % Make epochs
            
            % Calculate pepisode
            
            % Export data
            
        end % thisElec
        
    end % thisSession
    
end % thisSubject