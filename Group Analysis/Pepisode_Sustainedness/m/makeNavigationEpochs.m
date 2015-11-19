function makeNavigationEpochs(thisSubject, thisSession, thisElectrode, fixFile, findStoreFileName, thisEDF)
% This function will is used to create epochs from the time period when
% subjects are navigating in the correct arm of the maze, but have not yet
% entered the teleporter. These epochs will then be used to assess pepisode
% and compare them to pepisode during teleportation. The epochs are
% selected to be the earliest uninterrupted time period when the subject is
% in the correct arm of the maze. The epoch may not overlap with any
% portion of the teleportation epoch [-3 : +6 relative to teleporter
% entry].
%
% USAGE:    function makeNavigationEpochs(thisSubject, thisSession, thisElectrode, fixFile, findStoreFileName, thisEDF)
%
% INPUTS:
%   thisSubject: string for subjectID (e.g., 'UCDMC13')
%   thisSession: string for session (e.g., 'TeleporterA')
%   thisElectrode: string for the depth electrode of interest (e.g., 'LAD')
%   fixFile: 1 or 0 to indicate whether the findStore txt file from Unity
%       was manually fixed (this is only true for UCDMC15)
%   findStoreFileName: string for the filename of the unity txt file output
%       during navigation (e.g., 's1_patientTeleporterData.txt');
%
% OPTIONAL INPUTS (for UCDMC15 to account for 2 EDFs):
%   thisEDF: 1 or 2 to indicate which EDF for UCDMC15 (leave blank
%   otherwise)
%
% OUTPUTS:
%   This function will not return any outputs, but will save two files:
%   1. EEG dataset of the newly epoched data
%       'Epoched Data/' thisSubject '_' thisSession '_epoched_' thisElectrode '_navigation.set'
%   2. Mat file containing the list of real trial numbers for which
%   matching navigation epochs are available
%       'Mat Files/' thisSubject '_' thisSession '_' thisElectrode '_navigation_goodEpochs.mat'


%% TESTING
% clear all;
% thisSubject   = 'UCDMC15';
% thisSession   = 'TeleporterA';
% thisElectrode = 'LAD';
% fixFile       = 1;
% findStoreFileName = 's3_FindStore_TeleporterA_FIXED.txt';
% thisEDF = 1;



%% Load data

subjectDir = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' thisSubject '/'];

% load Mat file of good epochs
% var = goodEpochs
goodEpochsMatPath = [subjectDir 'Mat Files/' thisSubject '_' thisSession '_' thisElectrode '_noSpikes_noWaves_goodEpochs.mat'];
load(goodEpochsMatPath);

% load Mat file of epoch onset times
% var = epochsEDF
% N.B. epochsEDF1 & epochsEDF2 for UCDMC15
epochOnsetMatPath = [subjectDir 'Mat Files/' thisSubject '_' thisSession '_Epochs_Entry.mat'];
load(epochOnsetMatPath);

% load unepoched dataset with spikes and waves removed
% var = EEG.event.type (boundary vs. trial)
% var = EEG.event.latency
% var = EEG.event.duration
if ~exist('thisEDF', 'var')
    unepochedEEGPath = [subjectDir 'PreProcessing Intermediates/' thisSubject '_' thisSession '_unepoched_' thisElectrode '_noSpikes_noWaves.set'];
else
    unepochedEEGPath = [subjectDir 'PreProcessing Intermediates/' thisSubject '_' thisSession '_EDF' num2str(thisEDF) '_unepoched_' thisElectrode '_noSpikes_noWaves.set'];
end

EEG           = pop_loadset(unepochedEEGPath);
eventType     = {EEG.event.type};
boundaryLatency  = {EEG.event.latency};
eventDuration = {EEG.event.duration};

% load which trial numbers are contained within the unepoched EEG (some
% trials may have been removed if their onsets occurred within a boundary
% event)
unepochedEpochsPath = [subjectDir 'Mat Files/' thisSubject '_' thisSession '_' thisElectrode '_noSpikes_noWaves_unepoched_goodEpochs.mat'];
load(unepochedEpochsPath);

%% Clean data

% get timings in seconds of trial onset, entry into correct arm, and entry
% into the teleporter
[trialOnset, enterCorrectArm, enterTeleporter] = calculateTrialEventTiming([subjectDir 'Behavioral Data/' thisSession '/' findStoreFileName], fixFile, thisSession);

% restrict to good trials only
eventInd            = find(strcmpi('boundary', eventType) == 0);
enterTeleporterBins = cell2mat(boundaryLatency(eventInd));
boundaryInd         = find(strcmpi('boundary', eventType) == 1);
boundaryLatency     = cell2mat(boundaryLatency(boundaryInd));
boundaryDuration    = cell2mat(eventDuration(boundaryInd));

% goodNoWaveEpochs tells us which trials we still have in the unepoched EEG
% whereas goodEpochs tells us which trials we still have in the epoched EEG
% We'll now make all the data match what's in goodEpochs

% If one EDF (UCDMC13 & UCDMC14)
if ~exist('thisEDF', 'var')
    trialOnset      = trialOnset(goodEpochs);
    enterCorrectArm = enterCorrectArm(goodEpochs);
    enterTeleporter = enterTeleporter(goodEpochs);
    eTime           = eTime(goodEpochs);
    
    [~, ind, ~] = intersect(goodNoWaveEpochs, goodEpochs);
    enterTeleporterBins = enterTeleporterBins(ind);
else % UCDMC15
    if thisEDF == 1
        
        trialOnset(length(epochsEDF1) + 1:end)          = [];
        enterCorrectArm(length(epochsEDF1) + 1:end)     = [];
        enterTeleporter(length(epochsEDF1) + 1:end)     = [];
        eTime(length(epochsEDF1) + 1:end)               = [];
        goodEpochs(goodEpochs > length(epochsEDF1) + 1) = [];
        
        trialOnset      = trialOnset(goodEpochs);
        enterCorrectArm = enterCorrectArm(goodEpochs);
        enterTeleporter = enterTeleporter(goodEpochs);
        eTime           = eTime(goodEpochs);
        
        [~, ind, ~] = intersect(goodNoWaveEpochs1, goodEpochs);
        enterTeleporterBins = enterTeleporterBins(ind);
    else
        removeTrials = length(trialOnset) - length(epochsEDF2);
        removeInd    = find(goodEpochs <= removeTrials);
        
        goodEpochsEdf2 = goodEpochs;
        goodEpochsEdf2(removeInd) = [];
        
        trialOnset          = trialOnset(goodEpochsEdf2);
        enterCorrectArm     = enterCorrectArm(goodEpochsEdf2);
        enterTeleporter     = enterTeleporter(goodEpochsEdf2);
        eTime               = eTime(goodEpochsEdf2);
        
        
        
        goodEpochsReRef = goodEpochsEdf2 - removeTrials;
        [~, ind, ~] = intersect(goodNoWaveEpochs2, goodEpochsReRef);
        enterTeleporterBins = enterTeleporterBins(ind);
        
    end
end


% the teleportation epochs include a period before teleportation, which we
% do not want to overlap with, so calculate the # of bins for that here
preNtSec   = 1.830;
preFtSec   = 2.830;
preNtBins = round(preNtSec * EEG.srate);
preFtBins = round(preFtSec * EEG.srate);

% length of total epoch in bins (same as for normal teleporter epochs).
% We'll use this later to make sure we have enough contiguous timepoints.
epochBins = 9 * EEG.srate;
epochPreBins = 3 * EEG.srate; % how much to shift for time zero


% initialize output
epochTimeZero = nan(length(enterTeleporterBins), 1);
epochType     = epochTimeZero;
navTeleIntervalSecs = [];

for thisTrial = 1:length(enterTeleporterBins)
    
    % convert from seconds to bins
    thisEnterTeleSec = enterTeleporter(thisTrial);
    thisEnterTeleBin = enterTeleporterBins(thisTrial);
    thisEnterArmSec  = enterCorrectArm(thisTrial) - thisEnterTeleSec; % referenced to teleporter entry
    thisEnterArmBin  = round(thisEnterTeleBin + (thisEnterArmSec * EEG.srate));
    
    
    % identify the last possible valid time, which is the start of the
    % pre-teleportation part of the epoch
    if eTime(thisTrial) == 1 % NT
        thisEndBin = thisEnterTeleBin - preNtBins;
        epochType(thisTrial) = 1;
    else % FT
        thisEndBin = thisEnterTeleBin - preFtBins;
        epochType(thisTrial) = 2;
    end
    
    % check if there are any boundaries between thisEndBin and
    % thisEnterTeleBin. If so, update thisEndBin by accounting for the
    % missing time points
    thisBoundaryList = find(boundaryLatency > thisEndBin & boundaryLatency < thisEnterTeleBin);
    if ~isempty(thisBoundaryList)
        for thisBoundary = 1:length(thisBoundaryList)
            thisEndBin = thisEndArmBin + boundaryDuration(thisBoundaryList(thisBoundary));
        end
    end
    
    % check if there are any boundaries between thisEnterArmBin and
    % thisEnterTeleBin. If so, update thisEnterArmBin by accounting for the
    % missing time points
    thisBoundaryList = find(boundaryLatency > thisEnterArmBin & boundaryLatency < thisEnterTeleBin);
    if ~isempty(thisBoundaryList)
        for thisBoundary = 1:length(thisBoundaryList)
            thisEnterArmBin = thisEnterArmBin + boundaryDuration(thisBoundaryList(thisBoundary));
        end
    end
    
    
    
    % find the boundaries that occur between thisEnterArmBin and thisEndBin
    thisBoundaryList = find(boundaryLatency > thisEnterArmBin & boundaryLatency < thisEndBin);
    
    if length(thisBoundaryList) == 0 % no boundary events
        
        if thisEndBin - thisEnterArmBin >= epochBins
            epochTimeZero(thisTrial) = thisEnterArmBin + epochPreBins;
        else
            epochTimeZero(thisTrial) = nan;
        end
        
    else % there's at least 1 boundary
        
        % get the uninterrupted intervals
        thisBoundaryList = [thisEnterArmBin; boundaryLatency(thisBoundaryList)'; thisEndBin];
        theIntervals = diff(thisBoundaryList);
        
        % get intervals that are long enough
        goodIntervals = find(theIntervals >= epochBins);
        
        
        if length(goodIntervals) == 0 % if no intervals long enough
            epochTimeZero(thisTrial) = nan;
        else % pick the earliest one in time
            epochTimeZero(thisTrial) = thisBoundaryList(min(goodIntervals)) + epochPreBins;
        end
        
        
    end
    
    navTeleIntervalSecs(end + 1) = (thisEndBin - (epochTimeZero(thisTrial) + 2*epochPreBins))/EEG.srate;
    
end

% remove invalid trials
goodTrials    = find(~isnan(epochTimeZero));
epochTimeZero = epochTimeZero(goodTrials);
epochTimeZero = round(epochTimeZero);
epochType     = epochType(goodTrials);

if ~exist('thisEDF', 'var')
    goodNavEpochs = goodEpochs(goodTrials);
elseif thisEDF == 1
    goodNavEpochs = goodEpochs(goodTrials);
else
    goodNavEpochs = goodEpochsReRef(goodTrials) + removeTrials;
end

%% Epoch the data

% Insert events into EEG
thisEpoch = 1;
numSpikes = size(EEG.event,2);
for n = numSpikes+1:numSpikes+length(epochTimeZero)
    
    EEG.event(n).latency = epochTimeZero(thisEpoch);
    
    EEG.event(n).type = strcat('3', num2str(epochType(thisEpoch)));
    
    thisEpoch = thisEpoch + 1;
    
end


%% Save epoched data
if isempty(goodNavEpochs)
    sprintf(['\n\n\n\n No valid epochs found for ' thisSubject ' ' thisSession ' ' thisElectrode '.'])
else
    [newEEG, validEpochs] = pop_epoch(EEG,{'31', '32'},[-3 6]);
    goodNavEpochs = goodNavEpochs(validEpochs);
    if ~exist('thisEDF', 'var')
        pop_saveset(newEEG, [subjectDir 'Epoched Data/' thisSubject '_' thisSession '_epoched_' thisElectrode '_navigation.set']);
    else
        pop_saveset(newEEG, [subjectDir 'Epoched Data/' thisSubject '_' thisSession '_epoched_' thisElectrode '_EDF' num2str(thisEDF) '_navigation.set']);
    end
    
    % Save list of good epochs for this electrode
    % Also save list of intervals between end of navigation epoch and start
    % of teleportation epoch
    if ~exist('thisEDF', 'var')
        save([subjectDir 'Mat Files/' thisSubject '_' thisSession '_' thisElectrode '_navigation_goodEpochs.mat'],'goodNavEpochs');
        save([subjectDir 'Mat Files/' thisSubject '_' thisSession '_' thisElectrode '_navigation_teleportation_epoch_intervals.mat'], 'navTeleIntervalSecs');
    else
        save([subjectDir 'Mat Files/' thisSubject '_' thisSession '_' thisElectrode '_navigation_goodEpochs_EDF' num2str(thisEDF) '.mat'],'goodNavEpochs');
        save([subjectDir 'Mat Files/' thisSubject '_' thisSession '_' thisElectrode '_navigation_teleportation_epoch_intervals_EDF' num2str(thisEDF) '.mat'], 'navTeleIntervalSecs');
    end
        
end

