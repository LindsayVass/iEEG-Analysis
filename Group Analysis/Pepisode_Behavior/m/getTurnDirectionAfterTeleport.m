function getTurnDirectionAfterTeleport(inputPath, fixFile, teleporter, outputPath)
% This function will parse the raw text output from the Find Store part of
% the experiment and determine whether the patient turns the correct
% direction to face the target arm after exiting the teleporter. I am
% excluding trials for which the subject is already facing the correct arm
% when they exit the teleporter (i.e., the previous trial was for the store
% on the opposite arm). Note that trial numbers in the output are 1-indexed
% to match the EEG data.
%
% << getTurnDirectionAfterTeleport(inputPath, fixFile, teleporter, outputPath)
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
yRot        = txtdata{10};

%% extract outcome for each trial
% WrongArm = went down wrong arm of maze
% RightDir = turned the optimal direction
% WrongDir = turned the suboptimal direction
% Invalid = was already facing correct target

% initialize output
output = cell(1, 2);
output(1, :) = {'TrialNumber', 'Outcome'};

% get indices and times of each trial's onset
trialOnsetInds  = [1; find(diff(trialNumber)) + 1];

storeList    = {storeInfo.storeName};
oppStoreList = {storeInfo.oppositeStoreName};

for thisTrial = 1:length(trialOnsetInds)
    
    % no previous teleportation
    if thisTrial == 1
        output(end + 1, :) = {thisTrial, 'Invalid'};
        continue
    end
    
    % determine whether valid trial (previous store is not opposite store)
    thisStore    = target(trialOnsetInds(thisTrial) + 5); % add 5 because sometimes the first entry will mistakenly be the prev trial
    thisStoreInd = find(strcmpi(thisStore, storeList));
    
    if thisTrial == 2
        prevStore = target(1);
    else
        prevStore = target(trialOnsetInds(thisTrial - 1) + 5); % add 5 because sometimes the first entry will mistakenly be the prev trial
    end
    
    if strcmpi(prevStore, oppStoreList{thisStoreInd}) == 1
        % invalid trial
        output(size(output, 1) + 1, :) = {thisTrial, 'Invalid'};
        continue
    end
    
    % select the position data for this trial
    theTrialNumber = trialNumber(trialOnsetInds(thisTrial));
    trialInds      = find(trialNumber == theTrialNumber);
    trial_xPos     = xPos(trialInds);
    trial_zPos     = zPos(trialInds);
    
    % identify the first time point when patient crosses the correct boundary
    correctBoundaryCrossing = eval(['trial_' storeInfo(thisStoreInd).boundary]);
    correctBoundaryCross    = min(find(correctBoundaryCrossing));
    
    % identify the time points when patient crosses any of the incorrect
    % boundaries; if it's before the correct boundary crossing, then throw
    % out the trial
    wrongBoundInds = setdiff([1:1:4], thisStoreInd);
    wrongArm = 0;
    for thisBound = 1:length(wrongBoundInds)
        wrongBoundaryCrossing = eval(['trial_' storeInfo(wrongBoundInds(thisBound)).boundary]);
        wrongBoundaryCross    = min(find(wrongBoundaryCrossing));
        if wrongBoundaryCross < correctBoundaryCross
            output(end + 1, :) = {thisTrial, 'WrongArm'};
            wrongArm = 1;
            break
        end
    end
    if wrongArm == 1
        continue
    end
    
    %   output(end + 1, :) = {thisTrial, 'GoodTrial'};
    
    % Get yRot data while patient is still in central plaza
    centralInds = trialInds(1:correctBoundaryCross - 1);
    trial_yRot  = yRot(centralInds);
    
    % if valid trial, determine whether patient turned the correct
    % direction. Because of a coding error, we only have the y rotation in
    % quaternions not degrees. This means we can't actually distinguish
    % facing West vs East (though we can distinguish North vs South).
    tol = 0.05;
    switch thisStoreInd
        case 1 % Target to North; facing either East or West
            testDir = 1;
        case 2 % Target to East; facing either North or South
            northDiff = abs(trial_yRot(1));
            southDiff = abs(trial_yRot(1) - 1);
            if northDiff < southDiff % facing North
                testDir = 1;
            else
                testDir = 0;
            end
        case 3 % Target to South; facing either East or West
            testDir = 0;
        case 4 % Target to West; facing either North or South
            % Figure out current direction
            northDiff = abs(trial_yRot(1));
            southDiff = abs(trial_yRot(1) - 1);
            if northDiff < southDiff % facing North
                testDir = 1;
            else
                testDir = 0;
            end
    end
    
    wrongDir = find(abs(trial_yRot - testDir) < tol);
    if ~isempty(wrongDir)
        output(end + 1, :) = {thisTrial, 'WrongDir'};
    else
        output(end + 1, :) = {thisTrial, 'RightDir'};
    end
end

%% save the output
dlmcell(outputPath, output, 'delimiter', ',');