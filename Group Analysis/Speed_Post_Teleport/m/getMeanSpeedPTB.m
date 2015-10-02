% Get mean speed for each pre- and post-teleportation epoch

iEEGDir = '/Users/Lindsay/Documents/MATLAB/iEEG/';
load([iEEGDir 'Group Analysis/Subject Info/SessionInfo2.mat']);
saveFile = [iEEGDir 'Group Analysis/Speed_Post_Teleport/csv/UCDMC13_UCDMC14_speed_data.csv'];

%% Loop through subjects/sessions/electrodes
speedData = {'Subject', 'Session', 'Depth', 'TrialNumber', 'Interval', 'Speed'};
for thisSubject = 1:2
    
    for thisSession = 1:length(sessionInfo(thisSubject).teleporter)
        
        for thisDepth = 1:length(sessionInfo(thisSubject).teleporter(thisSession).depths)
            
            subjectID = sessionInfo(thisSubject).subjectID;
            sessionID = sessionInfo(thisSubject).teleporter(thisSession).name;
            depthID   = sessionInfo(thisSubject).teleporter(thisSession).depths(thisDepth).name;
            
            % Get final clean epoch numbers (only examine the trials we actually
            % analyzed)
            load(strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Mat Files/', subjectID{1}, '_', sessionID{1}, '_', depthID{1}, '_noSpikes_noWaves_goodEpochs.mat'))
            load(strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Mat Files/', subjectID{1}, '_', sessionID{1}, '_Epochs_Entry.mat'))
            
            epochsEDF = epochsEDF(goodEpochs);
            eTime     = eTime(goodEpochs);
            
            % Get interval times in ticks
            load(strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Mat Files/', subjectID{1}, '_', sessionID{1}, '_time_sync.mat'))
            epochsTicks = round((epochsEDF - time_sync_regression(2)) / time_sync_regression(1));
            
            
            % Load txt file containing position info
            if thisSession == 1
                fileName = strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Behavioral Data/', sessionID{1}, '/s', num2str(thisSubject), '_patientTeleporterData.txt');
            else
                fileName = strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Behavioral Data/', sessionID{1}, '/s', num2str(thisSubject), '_patientTeleporterData 2.txt');
            end
            fid = fopen(fileName);
            data =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
            fclose(fid);
            
            % Extract variables of interest
            systemTime = data{3};
            timeType   = data{7};
            xPos       = data{8};
            zPos       = data{9};
            
            % get sampling rate of the text file
            txtIntervalMs = round(mean(diff(systemTime)) / 10000); % 10000 ticks per millisecond
            
            
            for thisEpoch = 1:length(epochsTicks)
                
                % set interval length
                if eTime(thisEpoch) == 1
                    intLength = 1830;
                else
                    intLength = 2830;
                end
                
                teleEntry = epochsTicks(thisEpoch);
                teleExit  = teleEntry + 10000 * intLength;
                preStart  = teleEntry - 10000 * intLength;
                postEnd   = teleExit  + 10000 * intLength;
                
                % Extract position for each time point in each interval
                preStartInd  = findClosestTime(preStart, systemTime);
                teleEntryInd = findClosestTime(teleEntry, systemTime);
                teleExitInd  = findClosestTime(teleExit, systemTime);
                postEndInd   = findClosestTime(postEnd, systemTime);
                
                preInds  = preStartInd:teleEntryInd;
                postInds = teleExitInd:postEndInd;
                
                preX = xPos(preInds);
                preZ = zPos(preInds);
                
                postX = xPos(postInds);
                postZ = zPos(postInds);
                
                % Calculate mean speed for each epoch
                meanPre = calcSpeed(preX, preZ, txtIntervalMs);
                meanPost = calcSpeed(postX, postZ, txtIntervalMs);
                
                % Save data
                speedData(end+1, :) = {subjectID, sessionID, depthID, goodEpochs(thisEpoch), 'Pre', meanPre};
                speedData(end+1, :) = {subjectID, sessionID, depthID, goodEpochs(thisEpoch), 'Post', meanPost};
                
            end
        end
    end
end


dlmcell(saveFile, speedData, 'delimiter', ',')