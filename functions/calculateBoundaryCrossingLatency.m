function calculateBoundaryCrossingLatency(inputPath, fixFile, teleporter, outputPath)
% This function will parse the raw text output from the Find Store part of
% the experiment and calculate how long it takes the participant to enter
% the correct arm of the maze. Entry into the correct arm is defined as the
% time when they cross from the central plaza into the arm. I am excluding
% trials for which the subject is already facing the correct arm when they
% exit the teleporter (i.e., the previous trial was for the store on the
% opposite arm). Note that trial numbers in the output are 1-indexed to
% match the EEG data.
%
% USAGE: function calculateBoundaryCrossingLatency(inputPath, fixFile, teleporter, outputPath)
%
% INPUTS:
%   inputPath: path to the txt file output by Unity for the FindStore
%       portion of the experiment
%   fixFile: set to 1 if the txt file was manually fixed; 0 otherwise
%   teleporter: string, either 'TeleporterA' or 'TeleporterB'; selects the 
%       correct stores to use for the analysis
%   outputPath: path for the output csv file

%% parameters

% store details
storeInfo = struct('storeName', {}, 'oppositeStoreName', {}, 'boundary', {});

if strcmpi(teleporter, 'TeleporterA') == 1
    
    storeInfo(1).storeName         = 'Coffee Shop';
    storeInfo(1).oppositeStoreName = 'Grocery Store';
    storeInfo(1).boundary          = 'zPos > 586';
    
    storeInfo(2).storeName         = 'Florist';
    storeInfo(2).oppositeStoreName = 'Pet Store';
    storeInfo(2).boundary          = 'xPos > 536';
    
    storeInfo(3).storeName         = 'Grocery Store';
    storeInfo(3).oppositeStoreName = 'Coffee Shop';
    storeInfo(3).boundary          = 'zPos < 488';
    
    storeInfo(4).storeName         = 'Pet Store';
    storeInfo(4).oppositeStoreName = 'Florist';
    storeInfo(4).boundary          = 'xPos < 438';
    
else % TeleporterB
    storeInfo(1).storeName         = 'Ice Cream';
    storeInfo(1).oppositeStoreName = 'Clothing Store';
    storeInfo(1).boundary          = 'zPos > 586';
    
    storeInfo(2).storeName         = 'Music Store';
    storeInfo(2).oppositeStoreName = 'Dentist';
    storeInfo(2).boundary          = 'xPos > 536';
    
    storeInfo(3).storeName         = 'Clothing Store';
    storeInfo(3).oppositeStoreName = 'Ice Cream';
    storeInfo(3).boundary          = 'zPos < 488';
    
    storeInfo(4).storeName         = 'Dentist';
    storeInfo(4).oppositeStoreName = 'Music Store';
    storeInfo(4).boundary          = 'xPos < 438';
end



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
target      = txtdata{4};
xPos        = txtdata{8};
zPos        = txtdata{9};

%% extract latency for each trial

% initialize output
output = cell(1, 2);
output(1, :) = {'TrialNumber', 'Latency'};
output(2, :) = {1, NaN}; % first trial is undefined by definition

% get indices and times of each trial's onset
trialOnsetInds  = find(diff(trialNumber)) + 1;
trialOnsetTimes = timeElapsed(trialOnsetInds);

storeList    = {storeInfo.storeName};
oppStoreList = {storeInfo.oppositeStoreName};

for thisTrial = 1:length(trialOnsetInds)
    
    % determine whether valid trial (previous store is not opposite store)
    
    thisStore    = target(trialOnsetInds(thisTrial));
    thisStoreInd = find(strcmpi(thisStore, storeList));
    
    if thisTrial == 1
        prevStore = target(1);
    else
        prevStore = target(trialOnsetInds(thisTrial - 1));
    end
    
    if strcmpi(prevStore, oppStoreList{thisStoreInd}) == 1
        % invalid trial
        output(size(output, 1) + 1, :) = {thisTrial + 1, NaN};
        continue
    end
    
    % select the position data for this trial
    theTrialNumber = trialNumber(trialOnsetInds(thisTrial));
    trialInds      = find(trialNumber == theTrialNumber);
    trial_xPos     = xPos(trialInds);
    trial_zPos     = zPos(trialInds);
    
    % identify the first time point when patient crosses the boundary
    boundaryCrossings = eval(['trial_' storeInfo(thisStoreInd).boundary]);
    boundaryCross     = timeElapsed(trialInds(min(find(boundaryCrossings))));
    
    % calculate latency (boundaryCross - trialStart)
    trialStart = timeElapsed(trialOnsetInds(thisTrial));
    latency    = boundaryCross - trialStart;
    
    % output to summary
    output(size(output, 1) + 1, :) = {thisTrial + 1, latency};
    
end

%% save the output
dlmcell(outputPath, output, 'delimiter', ',');