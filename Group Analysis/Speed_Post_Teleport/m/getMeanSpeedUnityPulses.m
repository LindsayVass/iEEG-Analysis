% Get mean speed for each pre- and post-teleportation epoch

iEEGDir = '/Users/Lindsay/Documents/MATLAB/iEEG/';
load([iEEGDir 'Group Analysis/Subject Info/SessionInfo2.mat']);
saveFile = [iEEGDir 'Group Analysis/Speed_Post_Teleport/csv/UCDMC15_speed_data.csv'];

%% Loop through subjects/sessions/electrodes
speedData = {'Subject', 'Session', 'Depth', 'TrialNumber', 'Interval', 'Speed'};
for thisSubject = 3
    
    for thisSession = 1:length(sessionInfo(thisSubject).teleporter)
        
        for thisDepth = 1:length(sessionInfo(thisSubject).teleporter(thisSession).depths)
            
            subjectID = sessionInfo(thisSubject).subjectID;
            sessionID = sessionInfo(thisSubject).teleporter(thisSession).name;
            depthID   = sessionInfo(thisSubject).teleporter(thisSession).depths(thisDepth).name;
            
            % Get final clean epoch numbers (only examine the trials we actually
            % analyzed)
            load(strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Mat Files/', subjectID{1}, '_', sessionID{1}, '_', depthID{1}, '_noSpikes_noWaves_goodEpochs.mat'))
            load(strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Mat Files/', subjectID{1}, '_', sessionID{1}, '_Epochs_Entry.mat'))
            
            goodEpochs1 = goodEpochs(goodEpochs <= length(epochsEDF1));
            goodEpochs2 = goodEpochs(goodEpochs > length(epochsEDF1)) - (length(eType) - length(epochsEDF2));
            eTime       = eTime(goodEpochs);
            
            % Get interval times in ticks
            epochsTicks1 = ticksEDF1(goodEpochs1)';
            epochsTicks2 = ticksEDF2(goodEpochs2)';
            epochsTicks  = cat(1, epochsTicks1, epochsTicks2);
            
            
            % Load txt file containing position info
            fileName = strcat(iEEGDir, 'Subjects/', subjectID{1}, '/Behavioral Data/', sessionID{1}, '/s', num2str(thisSubject), '_FindStore_', sessionID{1}, '_FIXED.txt');
            fid = fopen(fileName);
            data =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n');
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