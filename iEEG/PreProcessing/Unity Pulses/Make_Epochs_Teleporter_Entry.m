% This script will load in the Unity behavioral output from the Patient
% Teleporter experiment and extract the timepoints of interest. The time
% point we're interested in is the time they enter the teleporter.

%% Set up script
clc;clear;close all;

addpath(genpath('/Users/Lindsay/MATLAB/functions/'));

% Paths
subjectDir    = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/';
unityDataPath = [subjectDir 'Behavioral Data/TeleporterB/s3_FindStore_TeleporterB_FIXED.txt'];
pulsesEDF1    = [subjectDir 'Raw Data/UCDMC15_TeleporterB_EDF1_pulse_timing.mat'];
pulsesEDF2    = [subjectDir 'Raw Data/UCDMC15_TeleporterB_EDF2_pulse_timing.mat'];
EEG_EDF1      = [subjectDir 'Raw Data/UCDMC15_TeleporterB_EDF1_badChansRemoved_trimmed.mat'];
EEG_EDF2      = [subjectDir 'Raw Data/UCDMC15_TeleporterB_EDF2_badChansRemoved_trimmed.mat'];

% save files
unepochedEEG1_savefile   = [subjectDir 'Raw Data/UCDMC15_TeleporterB_EDF1_unepoched.set'];
unepochedEEG2_savefile   = [subjectDir 'Raw Data/UCDMC15_TeleporterB_EDF2_unepoched.set'];
epochedEEG_savefile      = [subjectDir 'Epoched Data/UCDMC15_TeleporterB_epoched.set'];
epochs_saveFile          = [subjectDir 'Mat Files/UCDMC15_TeleporterB_Epochs_Entry.mat'];

% Epoch start/end times in seconds
eStart = -3;
eEnd = 6;

% bad chans -- these channels show spiking activity, so we'll remove them
% before re-referencing
badChans = {'RAD1' 'RAD2' 'RAD3' 'LAD1' 'LAD2'};


 %% Load data
 fid  = fopen(unityDataPath);
 data = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n'); % use this version for unity output that we fixed in Matlab
%  data =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1); % use this version for raw unity output

 fclose(fid);
 
 systemTime = data{3};
 target     = data{4};
 spaceType  = data{6};
 timeType   = data{7};
 xPos       = data{8};
 zPos       = data{9};
 yRot       = data{10};
 
 %% Identify time points of interest
 
 % first, we'll find the teleporter intervals 
 teleporterEntryInd = zeros(size(target));
 for i = 1:length(target) - 1
     
     % First, find the time point when the target changes from "Teleporter"
     % to a landmark
     t = strfind(target{i},'Teleporter');
     n = strfind(target{i+1},'Teleporter');
     if (isempty(t))
         continue
     elseif t > 0 && isempty(n) == 1
         
         % Get the xPos, zPos, and yRot from 5 samples before
         % teleportation event
         targetX = xPos(i-5);
         targetZ = zPos(i-5);
         targetY = yRot(i-5);
         
         % Find all samples between t0 and t-3 seconds that match that XYZ
         matchX = find(xPos == targetX);
         matchX(matchX < (i-150)) = [];
         
         matchZ = find(zPos == targetZ);
         matchZ(matchZ < (i-150)) = [];
         
         matchY = find(yRot == targetY);
         matchY(matchY < (i-150)) = [];
         
         % Find the earliest time point match
         minX = min(matchX);
         minZ = min(matchZ);
         minY = min(matchY);
         
         % Find the max among the minimum timepoints
         allMin = cat(1,minX,minZ,minY);
         trueEntry = max(allMin);
         
         teleporterEntryInd(trueEntry) = 1;
         
     end
     
     if i == length(target) - 1 % if it's the last teleportation event, there's no next target
         
         % Get the xPos, zPos, and yRot from 5 samples before
         % teleportation event
         targetX = xPos(end-5);
         targetZ = zPos(end-5);
         targetY = yRot(end-5);
         
         % Find all samples between t0 and t-3 seconds that match that XYZ
         matchX = find(xPos == targetX);
         matchX(matchX < (end-150)) = [];
         
         matchZ = find(zPos == targetZ);
         matchZ(matchZ < (end-150)) = [];
         
         matchY = find(yRot == targetY);
         matchY(matchY < (end-150)) = [];
         
         % Find the earliest time point match
         minX = min(matchX);
         minZ = min(matchZ);
         minY = min(matchY);
         
         % Find the max among the minimum timepoints
         allMin = cat(1,minX,minZ,minY);
         trueEntry = max(allMin);
         
         teleporterEntryInd(trueEntry) = 1;
         
     end
 end
 
 
 teleporterEntryInd = logical(teleporterEntryInd);
 
 %% Save teleporter entry epochs
 
 spaceTimeInd = teleporterEntryInd;
 spaceTimeInd = logical(spaceTimeInd);
 
 epochs = systemTime(spaceTimeInd);
 epochsSpace = spaceType(spaceTimeInd);
 epochsTime = timeType(spaceTimeInd);
 
 %% Convert epochs from ticks to EEG bins
 
 teleporterInds = find(spaceTimeInd);
 [~, epochsEDF1, epochsEDF2] = unityPulseTiming(systemTime(teleporterInds), pulsesEDF1, pulsesEDF2);
 
 %% For each epoch, set the trial type
 
 % SPACE
 % Near Space  = 1
 % Far Space   = 2
 
 % TIME
 % Short Time  = 1
 % Long Time   = 2
 
 % TYPE
 % Near Space Short Time = 1
 % Near Space Long Time  = 2
 % Far Space Short Time  = 3
 % Far Space Long Time   = 4
 
 eSpace = zeros(size(epochs));
 eTime  = eSpace;
 eType  = eSpace;
 
 for thisEpoch = 1:length(epochs)
     
     if strcmpi(epochsSpace(thisEpoch),'Near') == 1
         eSpace(thisEpoch) = 1;
     else
         eSpace(thisEpoch) = 2;
     end
     
     if strcmpi(epochsTime(thisEpoch),'Short') == 1
         eTime(thisEpoch) = 1;
     else
         eTime(thisEpoch) = 2;
     end
     
     if eSpace(thisEpoch) == 1 && eTime(thisEpoch) == 1
         eType(thisEpoch) = 1;
     elseif eSpace(thisEpoch) == 1 && eTime(thisEpoch) == 2
         eType(thisEpoch) = 2;
     elseif eSpace(thisEpoch) == 2 && eTime(thisEpoch) == 1
         eType(thisEpoch) = 3;
     elseif eSpace(thisEpoch) == 2 && eTime(thisEpoch) == 2
         eType(thisEpoch) = 4;
     end
         
         
 end
 
 % Save the results
 save(epochs_saveFile,'epochsEDF1','epochsEDF2','eSpace','eTime','eType');
 
 %% Re-reference the EEG data using all electrodes except marker
 
 % remove badChans
 load(EEG_EDF1);
 chanNames = {EEG.chanlocs.labels};
 ind = [];
 
 for thisChan = 1:length(badChans)
     ind(end+1) = find(strcmpi(chanNames,badChans(thisChan)) == 1);
     
 end
 
 EEG.data(ind,:) = [];
 EEG.nbchan = EEG.nbchan - length(ind);
 EEG.chanlocs(ind) = [];
 eeglab redraw;
 
 % re-reference
 refEEG = mean(EEG.data(1:end-1,:,:),1);
 refEEG = repmat(refEEG,[size(EEG.data,1)-1 1 1]);
 EEG.data(1:end-1,:,:) = EEG.data(1:end-1,:,:) - refEEG;
 
 % save the re-referenced data
 pop_saveset(EEG,unepochedEEG1_savefile);
 EEG1 = EEG;
 
 % do the same for EDF2
 % remove badChans
 load(EEG_EDF2);
 chanNames = {EEG.chanlocs.labels};
 ind = [];
 
 for thisChan = 1:length(badChans)
     ind(end+1) = find(strcmpi(chanNames,badChans(thisChan)) == 1);
     
 end
 
 EEG.data(ind,:) = [];
 EEG.nbchan = EEG.nbchan - length(ind);
 EEG.chanlocs(ind) = [];
 eeglab redraw;
 
 % re-reference
 refEEG = mean(EEG.data(1:end-1,:,:),1);
 refEEG = repmat(refEEG,[size(EEG.data,1)-1 1 1]);
 EEG.data(1:end-1,:,:) = EEG.data(1:end-1,:,:) - refEEG;
 
 % save the re-referenced data
 pop_saveset(EEG,unepochedEEG2_savefile);
 EEG2 = EEG;

%  %% Re-reference the EEG data separately for each strip/grid/depth
%  keyboard;
%  
%  % First we need to figure out which channels are on the same electrode
%  load(EEG_EDF1);
%  chanNames = char({EEG.chanlocs.labels});
%  chanNames = chanNames(1:end-1,1:3); % all except marker channel (last channel)
%  chanTypes = cellstr(unique(chanNames,'rows'));
%  chanNames = cellstr(chanNames);
%  
%  for thisType = 1:length(chanTypes)
%      
%      % find all channels on this electrode
%      ind = strcmpi(chanTypes(thisType),chanNames);
%      
%      % calculate mean EEG for all channels on this electrode and subtract
%      % it out
%      refEEG = mean(EEG.data(ind,:,:),1);
%      refEEG = repmat(refEEG,[sum(ind) 1 1]);
%      EEG.data(ind,:,:) = EEG.data(ind,:,:) - refEEG;
%      
%  end
%  
%  % save the re-referenced data
%  pop_saveset(EEG,unepochedEEG1_savefile);
%  EEG1 = EEG;
%  
%  % Now do the same thing for the second EDF file
%  load(EEG_EDF2);
%  chanNames = char({EEG.chanlocs.labels});
%  chanNames = chanNames(1:end-1,1:3); % all except marker channel (last channel)
%  chanTypes = cellstr(unique(chanNames,'rows'));
%  chanNames = cellstr(chanNames);
%  
%  for thisType = 1:length(chanTypes)
%      
%      % find all channels on this electrode
%      ind = strcmpi(chanTypes(thisType),chanNames);
%      
%      % calculate mean EEG for all channels on this electrode and subtract
%      % it out
%      refEEG = mean(EEG.data(ind,:,:),1);
%      refEEG = repmat(refEEG,[sum(ind) 1 1]);
%      EEG.data(ind,:,:) = EEG.data(ind,:,:) - refEEG;
%      
%  end
%  
%  % save the re-referenced data
%  pop_saveset(EEG,unepochedEEG2_savefile);
%  EEG2 = EEG;
 
%% Epoch the EEG data

% Insert events into EEG1
for n = 1:length(epochsEDF1)
    
    EEG1.event(n).latency = epochsEDF1(n);
    
    if eSpace(n) == 1 && eTime(n) == 1
        EEG1.event(n).type = '11';
    elseif eSpace(n) == 1 && eTime(n) == 2
        EEG1.event(n).type = '12';
    elseif eSpace(n) == 2 && eTime(n) == 1
        EEG1.event(n).type = '21';
    elseif eSpace(n) == 2 && eTime(n) == 2
        EEG1.event(n).type = '22';
    end
   
end

% Because of missing trials in the EEG recording break, identify which
% epoch is the first of the EDF2 file
eCount = length(eType) - length(epochsEDF2) + 1;

% Insert events into EEG2
for n = 1:length(epochsEDF2)
    
    EEG2.event(n).latency = epochsEDF2(n);
    
    if eSpace(eCount) == 1 && eTime(eCount) == 1
        EEG2.event(n).type = '11';
    elseif eSpace(eCount) == 1 && eTime(eCount) == 2
        EEG2.event(n).type = '12';
    elseif eSpace(eCount) == 2 && eTime(eCount) == 1
        EEG2.event(n).type = '21';
    elseif eSpace(eCount) == 2 && eTime(eCount) == 2
        EEG2.event(n).type = '22';
    end
   
    eCount = eCount + 1;
end

% Check if EEG data exists for the whole epoch of the trials that border
% the EDF file split
EEG1_eventEnd = EEG1.event(end).latency + EEG1.srate*eEnd;
if EEG1_eventEnd > size(EEG1.data,2)
    EEG1.event(end) = [];
end

EEG2_eventStart = EEG2.event(1).latency - EEG2.srate*eStart;
if EEG2_eventStart < 1
    EEG2.event(1) = [];
end

% Create epochs
[EEG1] = pop_epoch(EEG1,{},[eStart eEnd]);
[EEG2] = pop_epoch(EEG2,{},[eStart eEnd]);

%% reject bad epochs for each EEG

% manual rejection for EEG1
clear TMPREJ;
cmd = 'TMPREJ;';
eegplot(EEG1.data, 'winlength', 1, 'command', cmd, 'butlabel', 'Mark for Rejection','srate',EEG1.srate,'events',EEG1.event,'limits',[eStart*1000 eEnd*1000]);
readyNow = input('When finished rejecting epochs, press any key to continue to EEG2.','s');

if exist('TMPREJ', 'var')
    rejectedEEG1 = TMPREJ(:, 1);
    
    % convert to latency
    rejectedEEG1 = rejectedEEG1 - eStart*EEG1.srate + 1;
    
    % find the indices of the rejected trials
    rejectedEpochs1 = nan(size(rejectedEEG1));
    for thisReject = 1:length(rejectedEEG1)
        rejectedEpochs1(thisReject) = find(cell2mat({EEG1.event.latency}) == rejectedEEG1(thisReject));
    end
else
    rejectedEpochs1 = [];
end

% manual rejection for EEG2
clear TMPREJ;
eegplot(EEG2.data, 'winlength', 1, 'command', cmd, 'butlabel', 'Mark for Rejection','srate',EEG2.srate,'events',EEG2.event,'limits',[eStart*1000 eEnd*1000]);
readyNow = input('When finished rejecting epochs, press any key to continue.','s');

if exist('TMPREJ', 'var')
    rejectedEEG2 = TMPREJ(:, 1);
    
    % convert to latency
    rejectedEEG2 = rejectedEEG2 - eStart*EEG2.srate + 1;
    
    % find the indices of the rejected trials
    rejectedEpochs2 = nan(size(rejectedEEG2));
    for thisReject = 1:length(rejectedEEG2)
        rejectedEpochs2(thisReject) = find(cell2mat({EEG2.event.latency}) == rejectedEEG2(thisReject));
    end
else
    rejectedEpochs2 = [];
end

% save the indices of the rejected epochs
save(epochs_saveFile,'epochsEDF1','epochsEDF2','eSpace','eTime','eType','rejectedEpochs1','rejectedEpochs2');


%% merge the two data files
EEG = pop_mergeset(EEG1, EEG2);

% Save the epoched data
EEG = pop_saveset(EEG, 'filename', epochedEEG_savefile);
