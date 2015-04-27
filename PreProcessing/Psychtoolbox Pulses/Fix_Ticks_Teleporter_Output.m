% Fix tick timing for UCDMC14 Teleporter B. During the testing session, the
% laptop's clock updated the time, resulting in a 5 minute jump in the
% ticks recorded by Unity. This script will find that discontinuity,
% correct for it, and output a new version of the file with corrected
% ticks.
%
% Lindsay Vass 24 April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% set file path
fileName = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Behavioral Data/TeleporterB/s2_patientTeleporterData 2.txt';
saveFile = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Behavioral Data/TeleporterB/s2_patientTeleporterData_2_FIXED.txt';
saveMat = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Mat Files/teleporterB_tickjump.mat';


%% Parse the file
[trialNumber,timeElapsed,systemTime,target,armType,teleporterSpaceType,teleporterTimeType,xPos,zPos,yRot] = textread(fileName,'%n%n%n%s%s%s%s%n%n%n','delimiter',',','headerlines',1);

% save header for output
fid = fopen(fileName);
headerLine = fgets(fid);
fclose(fid);

% Calculate tick interval
tick_interval = diff(systemTime);

% Look at the data to see if there are any weird discontinuities
figure;
hist(tick_interval);
title('Interval between ticks');

% Find the time point at which the tick value jumps
jumpInd = find(tick_interval == max(tick_interval));

% Find the median time between ticks, not including the jump
intervalCopy = tick_interval;
intervalCopy(jumpInd) = [];
meanInterval = round(mean(intervalCopy));

% subtract out the jump from the system time
newSystemTime = systemTime;
jumpSubtract = tick_interval(jumpInd) - meanInterval;
newSystemTime(jumpInd:end) = newSystemTime(jumpInd:end) - jumpSubtract;
newInterval = diff(newSystemTime);

% Look at the new data
figure;
hist(newInterval);
title('NEW Interval between ticks');

%% Save the new file
fid = fopen(saveFile,'w');
fprintf(fid,'%s',headerLine);
for thisLine = 1:length(trialNumber)
    fprintf(fid,'%i,%f,%18.0f,%s,%s,%s,%s,%f,%f,%f\n',trialNumber(thisLine),timeElapsed(thisLine),newSystemTime(thisLine),target{thisLine},armType{thisLine},teleporterSpaceType{thisLine},teleporterTimeType{thisLine},xPos(thisLine),zPos(thisLine),yRot(thisLine));
end
fclose(fid);

% save the jump info so we can account for it in the psychtoolbox ticks
save(saveMat,'jumpSubtract');
