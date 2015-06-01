% UCDMC14 Preprocessing for Teleporter B.
% Lindsay Vass 27 April 2015

clear all; close all; clc;

%% Set up paths
subjectDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/';
subjectID  = 'UCDMC14';
teleporter  = 'TeleporterA';
edfFile    = 'UCDMC14_020415_teleporter.edf'; % in 'Raw Data' folder

% Psychtoolbox data
preTestPulses  = [subjectDir 'Raw Data/UCDMC14_020415_pre/SubjectUCDMC14_Data.mat'];
postTestPulses = [subjectDir 'Raw Data/UCDMC14_020415_post/SubjectUCDMC14_part2_Data_interim.mat'];

% Unity data
unityFindStorePath = [subjectDir 'Behavioral Data/' teleporter '/s2_patientTeleporterData.txt'];

% Epoch start/end times in seconds
eStart = -3;
eEnd = 6;

% we will pre-process the data separately for each depth electrode, so list
% them here
depthNames = {'LAD','LHD','RAD','RHD'};

% Two levels of "cleaning" the data. In the more liberal version
% (noSpikes), we remove dramatic perturbations in the EEG signal. In the
% more conservative version (noWaves), we also remove sharp waves.
labelsToRemove_noSpikes = {'spike','complex','other'};
labelsToRemove_noWaves = {'sharpWave'};

%% Set up filenames and paths

save_stem = [subjectID '_' teleporter];

addpath(genpath(subjectDir))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/iEEG/Amber Scripts/'))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b/'))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'))

% Make directories to put our data into
cd(subjectDir);
if ~exist([subjectDir 'PreProcessing Intermediates'],'dir')
    system('mkdir PreProcessing\ Intermediates');
end

if ~exist([subjectDir 'Epoched Data'],'dir')
    system('mkdir Epoched\ Data');
end

if ~exist([subjectDir 'Mat Files'],'dir')
    system('mkdir Mat\ Files');
end

if ~exist([subjectDir 'Figures'],'dir')
    system('mkdir Figures');
end

%% Load psychtoolbox data
load(preTestPulses);
pretest_master_event_list = master_event_list;

load(postTestPulses);
posttest_master_event_list = master_event_list;

%% Load patient LFP data in .edf format and put into EEG struct
cd ([subjectDir])

% Check if we've already done this stage of pre-processing
if exist(['PreProcessing Intermediates/' save_stem '_badChansRemoved.mat'],'file')
    % Ask user whether to use the old pulse timing file or re-calculate a
    % new one
    question = 'You already have a dataset with bad channels removed. What would you like to do?';
    questTitle = 'Redo pre-processing?';
    button1 = 'Use existing files';
    button2 = 'Redo pre-processing';
    answer = questdlg(question, questTitle, button1, button2, button2);
    
    if strcmpi(answer,button2) == 1
        redoBadChans = 1;
    else
        redoBadChans = 0;
    end
    
else
    redoBadChans = 1;
end

if redoBadChans == 1
    PutIntoEEG(subjectDir, subjectID, {edfFile}, 1, {save_stem});
end



%% Estimate time of pre- and post-test pulses and trim data before and after them

load(['PreProcessing Intermediates/' save_stem '_badChansRemoved.mat']);

% update marker channel number in case it changed
markerChanNumber = find(strcmpi('marker',{EEG.chanlocs.labels}));

% Check if we've already done this stage of pre-processing
if exist(['PreProcessing Intermediates/' save_stem '_badChansRemoved_trimmed.mat'],'file')
    % Ask user whether to use the old pulse timing file or re-calculate a
    % new one
    question = 'You already have a trimmed dataset. What would you like to do?';
    questTitle = 'Redo pre-processing?';
    button1 = 'Use existing files';
    button2 = 'Redo pre-processing';
    answer = questdlg(question, questTitle, button1, button2, button2);
    
    if strcmpi(answer, button2) == 1
        redoTrim = 1;
    else
        redoTrim = 0;
    end
    
else
    redoTrim = 1;
end

if redoTrim == 1
    % Plot the marker channel data and ask the user to identify the time of the
    % first pulse
    h = figure;
    x = [1:1:size(EEG.data,2)]/EEG.srate;
    plot(x,EEG.data(markerChanNumber,:));
    title('Marker channel activity');
    xlabel('Time (seconds)');
    prompt = '\n\nZoom in on the "Marker channel activity" plot until you can see the first pulse. \nInput the time in seconds just before the first pulse: ';
    firstPulseSec = input(prompt);
    firstPulseBin = firstPulseSec*EEG.srate;
    close(h);
    
    % Remove data before first pulse
    bufferSec = 5;
    bufferBin = bufferSec*EEG.srate;
    
    if firstPulseSec >= bufferSec % don't trim unless the pre-pulse EEG length exceeds the buffer
        maxTrim   = firstPulseBin - bufferBin - 1;
        trimmedBins = [1 maxTrim];
        EEG.data(:,1:maxTrim) = [];
        EEG.pnts = size(EEG.data,2);
        
        % Update time of the first pulse after trimming
        firstPulseBin = firstPulseBin - maxTrim - 1;
        firstPulseSec = bufferSec;
        
    end
    
    
    % Plot the marker channel data and ask the user to identify the time of the
    % first pulse of the post-test
    h = figure;
    x = [1:1:size(EEG.data,2)]/EEG.srate;
    plot(x,EEG.data(markerChanNumber,:));
    title('Marker channel activity');
    xlabel('Time (seconds)');
    prompt = '\n\nZoom in on the "Marker channel activity" plot until you can see the first pulse of the POST-TEST pulses. \nInput the time in seconds just before the first pulse: ';
    lastPulseSec = input(prompt);
    lastPulseBin = lastPulseSec*EEG.srate;
    close(h);
    
    
    eeglab redraw;
    
    % Save current state of EEG
    save(['PreProcessing Intermediates/' save_stem '_badChansRemoved_trimmed.mat'],'EEG','markerChanName','badChanNames','trimmedBins', 'firstPulseBin', 'lastPulseBin');
    
end


%% Calculate pulse timing synchronization
if exist([subjectDir 'Mat Files/' save_stem  '_time_sync.mat'],'file')
    % Ask user whether to use the old pulse timing file or re-calculate a
    % new one
    question = 'You already have pulse time synchronization for this dataset. What would you like to do?';
    questTitle = 'Redo pre-processing?';
    button1 = 'Use existing files';
    button2 = 'Redo pre-processing';
    answer = questdlg(question, questTitle, button1, button2, button2);
    
    if strcmpi(answer, button2) == 1
        redoPulses = 1;
    else
        redoPulses = 0;
    end
    
else
    redoPulses = 1;
end

if redoPulses == 1
    
    if ~exist('firstPulseBin','var') || ~exist('lastPulseBin','var')
        load(['PreProcessing Intermediates/' save_stem '_badChansRemoved_trimmed.mat']);
    end
    
    % Find all pulses
    prePulses = PulseFinderPTB(EEG, markerChanNumber, firstPulseBin, 'forward', pretest_master_event_list);
    postPulses = PulseFinderPTB(EEG, markerChanNumber, lastPulseBin, 'forward', posttest_master_event_list);
    
    % Get the timestamps from the stimulus computer for both sets of pulses
    pre_time_orig = pretest_master_event_list.time;
    post_time_orig = posttest_master_event_list.time;
    
    % Convert datenum time to ticks (this is what my code spits out from Unity,
    % plus it will be easier to fit a regression when the datetime is
    % represented as a single value).
    pre_time = nan(1,length(pre_time_orig));
    post_time = nan(1,length(post_time_orig));
    for i = 1:length(pre_time_orig)
        t = datenum(pre_time_orig{i});
        t = datenum2ticks(t);
        pre_time(i) = t;
    end
    
    for i = 1:length(post_time_orig)
        t = datenum(post_time_orig{i});
        t = datenum2ticks(t);
        post_time(i) = t;
    end
    
    % Calculate regression using both
    all_time = cat(2,pre_time,post_time);
    all_pulses = cat(1,prePulses,postPulses)';
    all_P = polyfit(all_time,all_pulses,1);
    all_y = polyval(all_P,all_time);
    
    % Plot the regression
    figure;
    plot(all_time,all_pulses,'k*')
    hold on;
    plot(all_time,all_y,'r-')
    title('All pulses')
    
    % Save regression values for later
    time_sync_regression = all_P;
    save([subjectDir 'Mat Files/' save_stem  '_time_sync.mat'],'time_sync_regression');
end

%% Find epochs of interest in the behavioral file

if exist([subjectDir 'Mat Files/' subjectID '_' teleporter '_Epochs_Entry.mat'],'file')
    % Ask user whether to use the old epoch timings or re-calculate new
    % ones
    question = 'You already have found the epoch timings for this dataset. What would you like to do?';
    questTitle = 'Redo pre-processing?';
    button1 = 'Use existing files';
    button2 = 'Redo pre-processing';
    answer = questdlg(question, questTitle, button1, button2, button2);
    
    if strcmpi(answer, button2) == 1
        redoFindEpochs = 1;
    else
        redoFindEpochs = 0;
    end
    
else
    redoFindEpochs = 1;
end

if redoFindEpochs == 1
    
    % Load in the time sync regression if it's not in the workspace
    if ~exist('time_sync_regression','var')
        load([subjectDir 'Mat Files/' save_stem  '_time_sync.mat'])
    end
    
    % Load in the unity output
    fid  = fopen(unityFindStorePath);
    data =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
    fclose(fid);
    
    % Extract variable of interest
    systemTime = data{3};
    target     = data{4};
    spaceType  = data{6};
    timeType   = data{7};
    xPos       = data{8};
    zPos       = data{9};
    yRot       = data{10};
    
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
    end % this Epoch
    
    % Convert to logicals
    teleporterEntryInd = logical(teleporterEntryInd);
    spaceTimeInd = teleporterEntryInd;
    spaceTimeInd = logical(spaceTimeInd);
    
    % Extract sytem time, space, and time types for each epoch
    epochs = systemTime(spaceTimeInd);
    epochsSpace = spaceType(spaceTimeInd);
    epochsTime = timeType(spaceTimeInd);
    
    % For each trial, convert time from ticks to bins
    teleporterInds = find(teleporterEntryInd);
    epochsEDF = nan(length(teleporterInds),1);
    
    % loop through trials
    for thisTrial = 1:length(teleporterInds)
        
        % get the time in ticks
        onsetTick = systemTime(teleporterInds(thisTrial));
        
        % convert to EEG bins
        onsetBin = onsetTick * time_sync_regression(1) + time_sync_regression(2);
        epochsEDF(thisTrial) = round(onsetBin);
        
    end % thisTrial
    
    % For each epoch, set the trial type
    
    %%%%%%SPACE %%%%%%%
    % Near Space  = 1 %
    % Far Space   = 2 %
    
    %%%%%% TIME %%%%%%%
    % Short Time  = 1 %
    % Long Time   = 2 %
    
    %%%%%%%%% TYPE %%%%%%%%%%%%%%
    % Near Space Short Time = 1 %
    % Near Space Long Time  = 2 %
    % Far Space Short Time  = 3 %
    % Far Space Long Time   = 4 %
    
    
    eSpace = zeros(size(epochsEDF));
    eTime  = eSpace;
    eType  = eSpace;
    
    for thisEpoch = 1:length(epochsEDF)
        
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
    save([subjectDir 'Mat Files/' subjectID '_' teleporter '_Epochs_Entry.mat'],'epochsEDF','eSpace','eTime','eType');
    
end

%% Re-reference the data using all electrodes except marker and spiking electrodes
if exist(['PreProcessing Intermediates/' save_stem '_unepoched.set'],'file')
    % Ask user whether to use the old re-referenced file or recalculate it
    question = 'You already have a re-referenced dataset. What would you like to do?';
    questTitle = 'Redo pre-processing?';
    button1 = 'Use existing files';
    button2 = 'Redo pre-processing';
    answer = questdlg(question, questTitle, button1, button2, button2);
    
    if strcmpi(button2,answer) == 1
        redoReRef = 1;
    else
        redoReRef = 0;
    end
else
    redoReRef = 1;
end

if redoReRef == 1
    
    load(['PreProcessing Intermediates/' save_stem '_badChansRemoved_trimmed.mat']);
    
    % select electrodes to use as reference
    refChans = 1:1:size(EEG.data,1);
    refChans(refChans == markerChanNumber) = [];
    for thisChan = 1:length(spikingChans)
        refChans(refChans == spikingChans(thisChan)) = [];
    end
    
    % re-reference
    refEEG = mean(EEG.data(refChans,:,:),1);
    refEEG = repmat(refEEG,[size(EEG.data,1)-1 1 1]);
    EEG.data(1:end-1,:,:) = EEG.data(1:end-1,:,:) - refEEG;
    
    % save the re-referenced data
    pop_saveset(EEG,['PreProcessing Intermediates/' save_stem '_unepoched.set']);
    
end

%% Separate the data for each depth electrode

% Check if we already have the data files
fileCheck = nan(length(depthNames),1);
for thisDepth = 1:length(depthNames)
    if exist(['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisDepth} '.set'],'file')
        fileCheck(thisDepth) = 1;
    else
        fileCheck(thisDepth) = 0;
    end
end

if isempty(find(fileCheck == 0))
    % Ask user whether to use the old re-referenced file or recalculate it
    question = 'You already have separate datasets for each depth electrode. What would you like to do?';
    questTitle = 'Redo pre-processing?';
    button1 = 'Use existing files';
    button2 = 'Redo pre-processing';
    answer = questdlg(question, questTitle, button1, button2, button2);
    
    if strcmpi(button2,answer) == 1
        redoDepth = 1;
    else
        redoDepth = 0;
    end
else
    redoDepth = 1;
end

if redoDepth == 1
    
    chanNames = char({EEG.chanlocs.labels});
    chanNames = chanNames(:,1:3);
    chanNames = cellstr(chanNames);
    
    EEG = pop_loadset(['PreProcessing Intermediates/' save_stem '_unepoched.set']);
    % Save a separate data set for each depth electrode
    for thisType = 1:length(depthNames)
        
        % find all channels on this electrode
        ind = strcmpi(depthNames(thisType),chanNames);
        
        % save a new dataset with just channels on this electrode
        EEGsub = EEG;
        EEGsub.data(ind==0,:,:) = [];
        EEGsub.chanlocs([ind == 0]) = [];
        EEGsub.nbchan = size(EEGsub.data,1);
        pop_saveset(EEGsub, ['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisType} '.set']);
        
    end
    
end


%% Annotate data with spikes and waves
% Load the data set for each electrode and mark time periods with spikes.
% Then set those time periods as NaN in the original dataset.

VISED_CONFIG = pop_edit_vised_config('vised_config.cfg');

for thisType = 1:length(depthNames)
    
    if exist(['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisType} '_noSpikes_noWaves.set'],'file')
        % Ask user whether to use the old pulse timing file or re-calculate a
        % new one
        question = 'You already have marked spikes for this channel. What would you like to do?';
        questTitle = 'Redo pre-processing?';
        button1 = 'Use existing files';
        button2 = 'Redo pre-processing';
        answer = questdlg(question, questTitle, button1, button2, button2);
        
        if strcmpi(answer,button1) == 1
            continue
        end
        
    else
        
        EEG = pop_loadset(['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisType} '.set']);
        EEG = eeg_checkset(EEG);
        eeglab redraw;
        
        % Initiate marks structure
        EEG.marks = marks_init(size(EEG.data), 0);
        EEG.marks = marks_add_label(EEG.marks, 'time_info', ...
            {'spike', [1, 0, 0], zeros(1, EEG.pnts)});
        EEG.marks = marks_add_label(EEG.marks, 'time_info', ...
            {'sharpWave', [0, 1, 0], zeros(1, EEG.pnts)});
        EEG.marks = marks_add_label(EEG.marks, 'time_info', ...
            {'complex', [0, 0, 1], zeros(1, EEG.pnts)});
        EEG.marks = marks_add_label(EEG.marks, 'time_info', ...
            {'other', [1, 0, 1], zeros(1, EEG.pnts)});
        
        % Annotate the data
        EEG = pop_vised(EEG, 'pop_gui', 'off');
        uiwait;
        
        % Save the marked dataset
        pop_saveset(EEG, ['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisType} '_marked.set']);

    end
end

%% Epoch the EEG data
load([subjectDir 'Mat Files/' subjectID '_' teleporter '_Epochs_Entry.mat']);
for thisDepth = 1:length(depthNames)
    
    % load the EEG data
    clear EEG;
    EEG = pop_loadset(['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisDepth} '_marked.set']);
    
    % Insert events into EEG
    thisEpoch = 1;
    numSpikes = size(EEG.event,2);
    for n = numSpikes+1:numSpikes+length(epochsEDF)
        
        EEG.event(n).latency = epochsEDF(thisEpoch);
        
        if eSpace(thisEpoch) == 1 && eTime(thisEpoch) == 1
            EEG.event(n).type = '11';
        elseif eSpace(thisEpoch) == 1 && eTime(thisEpoch) == 2
            EEG.event(n).type = '12';
        elseif eSpace(thisEpoch) == 2 && eTime(thisEpoch) == 1
            EEG.event(n).type = '21';
        elseif eSpace(thisEpoch) == 2 && eTime(thisEpoch) == 2
            EEG.event(n).type = '22';
        end
        
        thisEpoch = thisEpoch + 1;
        
    end
    
    % Save cleaned unepoched data
    % Save unepoched data with spikes removed
    indexes = marks_label2index(EEG.marks.time_info, labelsToRemove_noSpikes, 'indexes','exact','on');
    
    if isempty(indexes)
        warning('No spikes found. Saving the data as it is.')
        noSpikeEEG = pop_saveset(EEG, ['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisDepth} '_noSpikes.set']);
    else
        
        [noSpikeEEG,~] = pop_marks_select_data(EEG,'time marks',indexes);
        noSpikeEEG = pop_saveset(noSpikeEEG, ['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisDepth} '_noSpikes.set']);
    end
    
    % Save version with spikes and waves removed
    indexes = marks_label2index(noSpikeEEG.marks.time_info, labelsToRemove_noWaves, 'indexes', 'exact', 'on');
    
    if isempty(indexes)
        warning('No sharp waves found. Saving the data as it is.')
        noWaveEEG = pop_saveset(noSpikeEEG, ['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisDepth} '_noSpikes_noWaves.set']);
    else
        [noWaveEEG,~] = pop_marks_select_data(noSpikeEEG, 'time marks', indexes);
        noWaveEEG = pop_saveset(noWaveEEG, ['PreProcessing Intermediates/' save_stem '_unepoched_' depthNames{thisDepth} '_noSpikes_noWaves.set']);
    end
    
    %% Save cleaned epoched data with spikes removed
    [newEEG, goodEpochs] = pop_epoch(noSpikeEEG,{'11' '12' '21' '22'},[eStart eEnd]);
    
    % Save epoched data
    pop_saveset(newEEG, ['Epoched Data/' subjectID '_' teleporter '_epoched_' depthNames{thisDepth} '_noSpikes.set']);
    
    % Save list of good epochs for this electrode
    save(['Mat Files/' subjectID '_' teleporter '_' depthNames{thisDepth} '_noSpikes_goodEpochs.mat'],'goodEpochs');
    
    %% Save cleaned epoched data with spikes AND waves removed
    [newEEG, goodEpochs] = pop_epoch(noWaveEEG,{'11' '12' '21' '22'},[eStart eEnd]);
     
    % Save epoched data
    pop_saveset(newEEG, ['Epoched Data/' subjectID '_' teleporter '_epoched_' depthNames{thisDepth} '_noSpikes_noWaves.set']);
    
    % Save list of good epochs for this electrode
    save(['Mat Files/' subjectID '_' teleporter '_' depthNames{thisDepth} '_noSpikes_noWaves_goodEpochs.mat'],'goodEpochs');
end