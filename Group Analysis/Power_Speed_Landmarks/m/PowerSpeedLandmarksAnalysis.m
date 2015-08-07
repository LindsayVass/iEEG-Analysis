% The purpose of this script is to extract power, speed, and landmark
% richness values for each segment of navigation data. This script sets up
% the analyses, which are subsequently run by the two main functions:
%   runPowerSpeedLandmarksOneEDF
%   runPowerSpeedLandmarksTwoEDFs
%
% The only inputs that need to be provided are the length of the data
% segments in milliseconds (specified by "intervalMs") and the list of
% frequencies (specified by "frequencies").
%
% This script will output a comma separated text file for each
% subject/session/electrode, containing 7 columns:
%   Subject (e.g., UCDMC13)
%   Session (e.g., TeleporterA)
%   Electrode (e.g., LAD1)
%   Frequency (in Hz)
%   Power (mean of the log(power) for each segment of data)
%   Speed (mean speed for each segment of data; has been prewhitened)
%   Landmarks (3 conditions: rich arm (many landmarks), poor arm (few
%       landmarks), or central plaza)
%
% Author: Lindsay Vass
% Date: 7 August 2015


% parameters
intervalMs  = 200; % length of the data segments in milliseconds
frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % log-spaced frequencies to use for analysis
subjectsDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';

load('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo.mat');

% behavioral files are all named differently, so add that info here
sessionInfo(1).behavioral = {'s1_patientTeleporterData.txt'};
sessionInfo(2).behavioral = {'s2_patientTeleporterData.txt', 's2_patientTeleporterData 2.txt'};
sessionInfo(3).behavioral = {'s3_FindStore_TeleporterA_FIXED.txt', 's3_FindStore_TeleporterB_FIXED.txt'};

for thisSubject = 3:size(sessionInfo, 2)
    
    subjectID = sessionInfo(thisSubject).subjectID;
    
    for thisSession = 1:size(sessionInfo(thisSubject).teleporter, 2)
        
        teleporter = sessionInfo(thisSubject).teleporter{thisSession};
        
        for thisElectrode = 1:size(sessionInfo(thisSubject).chanList, 2)
            
            electrode  = sessionInfo(thisSubject).chanList{thisElectrode};
            numEDFs    = sessionInfo(thisSubject).numEDFs(thisSession);
            behavioral = sessionInfo(thisSubject).behavioral{thisSession}; 
            
            if numEDFs == 1
                behavioralPath   = [subjectsDir subjectID '/Behavioral Data/' teleporter '/' behavioral];
                timeSyncPath     = [subjectsDir subjectID '/Mat Files/' subjectID '_' teleporter '_time_sync.mat'];
                unepochedEEGPath = [subjectsDir subjectID '/PreProcessing Intermediates/' subjectID '_' teleporter '_unepoched_' electrode(1:3) '_marked.set'];
                [meanLogPowerVals, speedWhite, landmarkClean] = runPowerSpeedLandmarksOneEDF(behavioralPath, timeSyncPath, unepochedEEGPath, frequencies, electrode, intervalMs);
            else
                behavioralPath    = [subjectsDir subjectID '/Behavioral Data/' teleporter '/' behavioral];
                pulseTimingPaths  = {[subjectsDir subjectID '/Mat Files/' subjectID '_' teleporter '_EDF1_pulse_timing.mat'], ...
                                     [subjectsDir subjectID '/Mat Files/' subjectID '_' teleporter '_EDF2_pulse_timing.mat']};
                unepochedEEGPaths = {[subjectsDir subjectID '/PreProcessing Intermediates/' subjectID '_' teleporter '_EDF1_unepoched_' electrode(1:3) '_marked.set'], ...
                                     [subjectsDir subjectID '/PreProcessing Intermediates/' subjectID '_' teleporter '_EDF2_unepoched_' electrode(1:3) '_marked.set']};
                [meanLogPowerVals, speedWhite, landmarkClean] = runPowerSpeedLandmarksTwoEDFs(behavioralPath, pulseTimingPaths, unepochedEEGPaths, frequencies, electrode, intervalMs);
            end
            
            powerVector     = reshape(meanLogPowerVals', [], 1);
            frequencyVector = repmat(frequencies, [size(meanLogPowerVals, 2) 1]);
            frequencyVector = reshape(frequencyVector, [], 1);
            speedVector     = repmat(speedWhite, [length(frequencies) 1]);
            landmarkVector  = repmat(landmarkClean, [length(frequencies) 1]);
            landmarkVector  = [landmarkVector{:}]';
            subjectVector   = repmat(subjectID, size(powerVector));
            sessionVector   = repmat(teleporter, size(powerVector));
            electrodeVector = repmat(electrode, size(powerVector));
            
            thisTable = table(subjectVector, sessionVector, electrodeVector, frequencyVector, powerVector, speedVector, landmarkVector);
            thisTable.Properties.VariableNames = {'Subject', 'Session', 'Electrode', 'Frequency', 'Power', 'Speed', 'Landmarks'};
            
            savePath = ['/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power_Speed_Landmarks/mat/' subjectID '_' teleporter '_' electrode '.txt'];
            writetable(thisTable, savePath);
            
        end
    end
end
