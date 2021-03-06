
clear all;
addpath('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Sustainedness/m/');
load('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo.mat');

%% update to add file names

% UCDMC13
sessionInfo(1).findStoreFileName = {'s1_patientTeleporterData.txt'};
sessionInfo(1).Electrode = {'LAD', 'LHD'};

% UCDMC14
sessionInfo(2).findStoreFileName = {'s2_patientTeleporterData.txt', 's2_patientTeleporterData 2.txt'};
sessionInfo(2).Electrode = {'LAD', 'LHD', 'RAD', 'RHD'};

% UCDMC15
sessionInfo(3).findStoreFileName = {'s3_FindStore_TeleporterA_FIXED.txt', 's3_FindStore_TeleporterB_FIXED.txt'};
sessionInfo(3).Electrode = {'LAD', 'LHD', 'RAD', 'RHD'};

%% run function

for thisSubject = 1:size(sessionInfo, 2)
    
    subjectID = sessionInfo(thisSubject).subjectID;
    
    for thisSession = 1:size(sessionInfo(thisSubject).teleporter, 2)
        
        sessionID = sessionInfo(thisSubject).teleporter{thisSession};
        findStoreFileName = sessionInfo(thisSubject).findStoreFileName{thisSession};
        
        if ~isempty(strfind(sessionInfo(thisSubject).findStoreFileName(thisSession), 'FIXED'))
            fixFile = 1;
        else
            fixFile = 0;
        end
        
        for thisElectrode = 1:size(sessionInfo(thisSubject).Electrode, 2)
            
            electrodeID = sessionInfo(thisSubject).Electrode{thisElectrode};
            
            for thisEDF = 1:sessionInfo(thisSubject).numEDFs
                
                if sessionInfo(thisSubject).numEDFs > 1
                    
                    makeNavigationEpochs(subjectID, sessionID, electrodeID, fixFile, findStoreFileName, thisEDF);
                    
                else
                    
                    makeNavigationEpochs(subjectID, sessionID, electrodeID, fixFile, findStoreFileName);
                    
                end
                
            end
            
            if sessionInfo(thisSubject).numEDFs > 1
                
                EEG1path = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/Epoched Data/' subjectID '_' sessionID '_epoched_' electrodeID '_EDF1_navigation.set'];
                EEG2path = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/Epoched Data/' subjectID '_' sessionID '_epoched_' electrodeID '_EDF2_navigation.set'];
                saveEEG = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/Epoched Data/' subjectID '_' sessionID '_epoched_' electrodeID '_navigation.set'];
                
                if exist(EEG1path, 'file') == 0 & exist(EEG2path, 'file') ~= 0
                    EEG2 = pop_loadset(EEG2path);
                    pop_saveset(EEG2, saveEEG);
                elseif exist(EEG1path, 'file') ~= 0 & exist(EEG2path, 'file') == 0
                    EEG1 = pop_loadset(EEG1path);
                    pop_saveset(EEG1, saveEEG);
                elseif exist(EEG1path, 'file') ~=0 & exist(EEG2path, 'file') ~= 0
                    EEG1 = pop_loadset(EEG1path);
                    EEG2 = pop_loadset(EEG2path);
                    mergeEEG = pop_mergeset(EEG1, EEG2);
                    pop_saveset(mergeEEG, saveEEG);
                else
                    warning(['No EEG data found for ' subjectID ' ' sessionID ' ' electrodeID '.'])
                end
                
            end
            
        end
        
        
    end
end