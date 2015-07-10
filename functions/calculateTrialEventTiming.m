function [trialOnset, enterCorrectArm, enterTeleporter] = calculateTrialEventTiming(inputPath, fixFile, teleporter)
% This function will parse the raw text output from the Find Store part of
% the experiment and return three times:
%   1. Trial onset
%   2. Entry into correct arm of the maze
%   3. Entry into teleporter
%
% USAGE: [trialOnset, enterCorrectArm, enterTeleporter] = calculateTrialEventTiming(inputPath, fixFile, teleporter)
%
% INPUTS:
%   inputPath: path to the txt file output by Unity for the FindStore
%       portion of the experiment
%   fixFile: set to 1 if the txt file was manually fixed; 0 otherwise
%   teleporter: string, either 'TeleporterA' or 'TeleporterB'; selects the
%       correct stores to use for the analysis
%
% OUTPUT:
%   trialOnset: time in seconds of trial onset (where 0 is beginning of
%       experiment)
%   enterCorrectArm: time in seconds at which the subject entered the
%       correct arm of the maze (0 is beginning of experiment)
%   enterTeleporter: time in seconds at which the subject entered the
%       teleporter (0 is beginning of experiment)

%% parameters

% store details
storeInfo = struct('storeName', {}, 'boundary', {});

if strcmpi(teleporter, 'TeleporterA') == 1
    
    storeInfo(1).storeName         = 'Coffee Shop';
    storeInfo(1).boundary          = 'zPos > 586';
    
    storeInfo(2).storeName         = 'Florist';
    storeInfo(2).boundary          = 'xPos > 536';
    
    storeInfo(3).storeName         = 'Grocery Store';
    storeInfo(3).boundary          = 'zPos < 488';
    
    storeInfo(4).storeName         = 'Pet Store';
    storeInfo(4).boundary          = 'xPos < 438';
    
else % TeleporterB
    storeInfo(1).storeName         = 'Ice Cream';
    storeInfo(1).boundary          = 'zPos > 586';
    
    storeInfo(2).storeName         = 'Music Store';
    storeInfo(2).boundary          = 'xPos > 536';
    
    storeInfo(3).storeName         = 'Clothing Store';
    storeInfo(3).boundary          = 'zPos < 488';
    
    storeInfo(4).storeName         = 'Dentist';
    storeInfo(4).boundary          = 'xPos < 438';
end

storeList    = {storeInfo.storeName};

%% parse input

% load the unity text file
fid  = fopen(inputPath);

if fixFile == 1
    txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n');
else
    txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
end

fclose(fid);

trialNumber = txtdata{1};
timeElapsed = txtdata{2};
timeTicks   = txtdata{3};
target      = txtdata{4};
xPos        = txtdata{8};
zPos        = txtdata{9};
yRot        = txtdata{10};

%% extract timings for each trial

% initialize output
trialOnset      = nan(64, 1);
enterCorrectArm = nan(64, 1);
enterTeleporter = nan(64, 1);

% get times of each trial's onset
trialOnsetInds  = [1; find(diff(trialNumber)) + 1];
trialOnset      = timeElapsed(trialOnsetInds);

% get times of correct arm entry
for thisTrial = 1:length(trialOnsetInds)
    
    thisStore    = target(trialOnsetInds(thisTrial) + 5); % sometimes it registers an incorrect target on the first time point
    thisStoreInd = find(strcmpi(thisStore, storeList));
    
    % select the position data for this trial
    theTrialNumber = trialNumber(trialOnsetInds(thisTrial));
    trialInds      = find(trialNumber == theTrialNumber);
    trial_xPos     = xPos(trialInds);
    trial_zPos     = zPos(trialInds);
    
    % identify the first time point when patient crosses the boundary
    boundaryCrossings = eval(['trial_' storeInfo(thisStoreInd).boundary]);
    boundaryCross     = timeElapsed(trialInds(min(find(boundaryCrossings))));
    
    try
        enterCorrectArm(thisTrial) = boundaryCross;
        
    catch % in case subject did not complete this trial
        enterCorrectArm(thisTrial) = nan;
    end
    
    
end

% get times of teleporter entry
thisTrial = 1;
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
        
        enterTeleporter(thisTrial) = timeElapsed(trueEntry);
        thisTrial = thisTrial + 1;
        
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
        
        enterTeleporter(thisTrial) = timeElapsed(trueEntry);
        
    end
end 