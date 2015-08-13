% This script will perform the full pre-processing for EEG data obtained
% with unity pulses and two EDF files
%
% Lindsay Vass 14 May 2015

clear all;close all;clc;

%% Define parameters for pre-processing
subject_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/';
subject_id  = 'UCDMC15';
teleporter  = 'TeleporterA';
edf_file_1  = [subject_dir 'Raw Data/UCDMC15_04_01_15_teleporter.edf'];
edf_file_2  = [subject_dir 'Raw Data/UCDMC15_04_01_15_teleporter1.edf'];

% txt files output by Unity
unityPulsesFile = [subject_dir 'Behavioral Data/TeleporterA/s3_PULSES_.txt'];
unityDataPath   = [subject_dir 'Behavioral Data/TeleporterA/s3_FindStore_TeleporterA_FIXED.txt'];
fixFile         = 1; % set to 1 if the FindStore txt file was manually fixed; set to 0 otherwise

% we will pre-process the data separately for each depth electrode, so list
% them here
depthNames = {'LAD','LHD','RAD','RHD'};

% Epoch start/end times in seconds
eStart = -3;
eEnd = 6;

labelsToRemove_noSpikes = {'spike','complex','other'};
labelsToRemove_noWaves = {'spike','complex','other','sharpWave'};
%% Set up filenames and paths
edf_files = {edf_file_1, edf_file_2};

save_stem_1 = [subject_id '_' teleporter '_EDF1'];
save_stem_2 = [subject_id '_' teleporter '_EDF2'];
save_stems  = {save_stem_1, save_stem_2};

pulseSaveFile = [subject_dir 'Mat Files/' subject_id '_' teleporter '_pulses.mat'];

addpath(genpath(subject_dir))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/iEEG/Amber Scripts/'))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b/'))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'))

% Make directories to put our data into
cd(subject_dir);
if ~exist([subject_dir 'PreProcessing Intermediates'],'dir')
    system('mkdir PreProcessing\ Intermediates');
end

if ~exist([subject_dir 'Epoched Data'],'dir')
    system('mkdir Epoched\ Data');
end

if ~exist([subject_dir 'Mat Files'],'dir')
    system('mkdir Mat\ Files');
end

if ~exist([subject_dir 'Figures'],'dir')
    system('mkdir Figures');
end


%% Load EDF data and put into EEGLAB format
for thisEDF = 1:length(edf_files)
    
    % Check if we've already done this stage of pre-processing
    if exist(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved.mat'],'file')
        % Ask user whether to use the old pulse timing file or re-calculate a
        % new one
        question = 'You already have a dataset with bad channels removed. What would you like to do?';
        questTitle = 'Redo pre-processing?';
        button1 = 'Use existing files';
        button2 = 'Redo pre-processing';
        answer = questdlg(question, questTitle, button1, button2, button2);
        
        if strcmpi(answer,button2) == 1
            PutIntoEEG(subject_dir, subject_id, edf_files, thisEDF, save_stems);
        end
        
    else
        if thisEDF == 1
            PutIntoEEG(subject_dir, subject_id, edf_files, thisEDF, save_stems);
        else
            PutIntoEEG(subject_dir, subject_id, edf_files, thisEDF, save_stems);
        end
    end
end

%% Find pulses
for thisEDF = 1:length(edf_files)
    
    load(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved.mat']);
    
    if exist(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'file')
        
        % Ask user whether to use the old pulse timing file or re-calculate a
        % new one
        question = 'You already have a trimmed dataset. What would you like to do?';
        questTitle = 'Redo pre-processing?';
        button1 = 'Use existing files';
        button2 = 'Redo pre-processing';
        answer = questdlg(question, questTitle, button1, button2, button2);
        
        if strcmpi(answer,button2) == 1
            
            if thisEDF == 1
                %% Find the first pulse of EDF1
                
                % update marker channel number in case it changed
                markerChanNumber = find(strcmpi('marker',{EEG.chanlocs.labels}));
                
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
                maxTrim   = firstPulseBin - bufferBin - 1;
                trimmedBins = [1 maxTrim];
                
                EEG.data(:,1:maxTrim) = [];
                EEG.pnts = size(EEG.data,2);
                
                % Save current state of EEG
                save(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'EEG','markerChanName','badChanNames','trimmedBins');
                
            else
                %% Find the last pulse of EDF2
                
                % update marker channel number in case it changed
                markerChanNumber = find(strcmpi('marker',{EEG.chanlocs.labels}));
                
                % Plot the marker channel data and ask the user to identify the time of the
                % last pulse
                h = figure;
                x = [1:1:size(EEG.data,2)]/EEG.srate;
                plot(x,EEG.data(markerChanNumber,:));
                title('Marker channel activity');
                xlabel('Time (seconds)');
                prompt = '\n\nZoom in on the "Marker channel activity" plot until you can see the LAST pulse. \nInput the time in seconds just after the last pulse: ';
                lastPulseSec = input(prompt);
                lastPulseBin = lastPulseSec*EEG.srate;
                close(h);
                
                % Remove data after last pulse
                bufferSec = 5;
                bufferBin = bufferSec*EEG.srate;
                maxTrim   = lastPulseBin + bufferBin + 1;
                trimmedBins = [maxTrim size(EEG.data,2)];
                
                EEG.data(:,maxTrim:end) = [];
                EEG.pnts = size(EEG.data,2);
                
                % Save current state of EEG
                save(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'EEG','markerChanName','badChanNames','trimmedBins');
            end
        end
        
    else
        if thisEDF == 1
            %% Find the first pulse of EDF1
            
            % update marker channel number in case it changed
            markerChanNumber = find(strcmpi('marker',{EEG.chanlocs.labels}));
            
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
            maxTrim   = firstPulseBin - bufferBin - 1;
            trimmedBins = [1 maxTrim];
            
            EEG.data(:,1:maxTrim) = [];
            EEG.pnts = size(EEG.data,2);
            
            % Save current state of EEG
            save(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'EEG','markerChanName','badChanNames','trimmedBins','maxTrim','firstPulseBin');
            
        else
            %% Find the last pulse of EDF2
            
            % update marker channel number in case it changed
            markerChanNumber = find(strcmpi('marker',{EEG.chanlocs.labels}));
            
            % Plot the marker channel data and ask the user to identify the time of the
            % last pulse
            h = figure;
            x = [1:1:size(EEG.data,2)]/EEG.srate;
            plot(x,EEG.data(markerChanNumber,:));
            title('Marker channel activity');
            xlabel('Time (seconds)');
            prompt = '\n\nZoom in on the "Marker channel activity" plot until you can see the LAST pulse. \nInput the time in seconds just after the last pulse: ';
            lastPulseSec = input(prompt);
            lastPulseBin = lastPulseSec*EEG.srate;
            close(h);
            
            % Remove data after last pulse
            bufferSec = 5;
            bufferBin = bufferSec*EEG.srate;
            maxTrim   = lastPulseBin + bufferBin + 1;
            trimmedBins = [maxTrim size(EEG.data,2)];
            
            EEG.data(:,maxTrim:end) = [];
            EEG.pnts = size(EEG.data,2);
            
            % Save current state of EEG
            save(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'EEG','markerChanName','badChanNames','trimmedBins','maxTrim','lastPulseBin');
        end
    end
    
    
    
    %% Find all pulses
    if exist(['Mat Files/' save_stems{thisEDF} '_pulse_timing.mat'],'file')
        
        % Ask user whether to use the old pulse timing file or re-calculate a
        % new one
        question = 'A pulse timing file already exists for this subject.';
        questTitle = 'Re-calculate pulse timing?';
        button1 = 'Use existing pulses';
        button2 = 'Re-calculate pulses';
        answer = questdlg(question, questTitle, button1, button2, button2);
        
        if strcmpi(answer, button1) == 1
            
        elseif strcmpi(answer, button2) == 1
            
            if thisEDF == 1
                
                if ~exist('firstPulseBin','var')
                    load(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'firstPulseBin','maxTrim');
                end
                
                % Update time of the first pulse
                firstPulseBin = firstPulseBin - maxTrim - 1;
                firstPulseSec = bufferSec;
                
                PulseFinder(subject_dir, subject_id, teleporter, EEG, markerChanNumber, firstPulseBin, 'forward', 250);
                
            else % thisEDF == 2
                
                
                if ~exist('lastPulseBin','var')
                    load(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'lastPulseBin');
                end
                
                PulseFinder(subject_dir, subject_id, teleporter, EEG, markerChanNumber, lastPulseBin, 'backward', 250);
                
            end
        end
        
    else
        if thisEDF == 1
            
            if ~exist('firstPulseBin','var')
                load(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'firstPulseBin','maxTrim');
            end
            
            % Update time of the first pulse
            firstPulseBin = firstPulseBin - maxTrim - 1;
            firstPulseSec = bufferSec;
            
            PulseFinder(subject_dir, subject_id, teleporter, EEG, markerChanNumber, firstPulseBin, 'forward', 250);
            
        else % thisEDF == 2
            
            if ~exist('lastPulseBin','var')
                load(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat'],'lastPulseBin');
            end
            
            PulseFinder(subject_dir, subject_id, teleporter, EEG, markerChanNumber, lastPulseBin, 'backward', 250);
            
        end
    end
end

%% Re-reference the data using all electrodes except marker and spiking electrodes
for thisEDF = 1:length(edf_files)
    
    if exist(['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched.set'],'file')
        % Ask user whether to use the old pulse timing file or re-calculate a
        % new one
        question = 'You already have a re-referenced dataset. What would you like to do?';
        questTitle = 'Redo pre-processing?';
        button1 = 'Use existing files';
        button2 = 'Redo pre-processing';
        answer = questdlg(question, questTitle, button1, button2, button2);
        
        if strcmpi(button1,answer) == 1
            continue
        end
    else
        load(['PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved_trimmed.mat']);
        
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
        pop_saveset(EEG,['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched.set']);
        
    end
    
end % thisEDF

%% Separate the data for each depth electrode
chanNames = char({EEG.chanlocs.labels});
chanNames = chanNames(:,1:3);
chanNames = cellstr(chanNames);


for thisEDF = 1:length(edf_files)
    EEG = pop_loadset(['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched.set']);
    % Save a separate data set for each depth electrode
    for thisType = 1:length(depthNames)
        
        % find all channels on this electrode
        ind = strcmpi(depthNames(thisType),chanNames);
        
        % save a new dataset with just channels on this electrode
        EEGsub = EEG;
        EEGsub.data(ind==0,:,:) = [];
        EEGsub.chanlocs([ind == 0]) = [];
        EEGsub.nbchan = size(EEGsub.data,1);
        pop_saveset(EEGsub, ['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisType} '.set']);
        
    end
    
end % thisEDF

%% Annotate data with spikes and waves
for thisEDF = 1:length(edf_files)
    
    % Load the data set for each electrode and mark time periods with spikes.
    % Then set those time periods as NaN in the original dataset.
    if thisEDF == 1
        VISED_CONFIG = pop_edit_vised_config('vised_config.cfg');
    end
    for thisType = 1:length(depthNames)
        
        if exist(['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisType} '_noSpikes_noWaves.set'],'file')
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
            
            EEG = pop_loadset(['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisType} '.set']);
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
            pop_saveset(EEG, ['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisType} '_marked.set']);
            
        end
    end
    
end % thisEDF



%% EPOCH the data

if exist([subject_dir 'Mat Files/' subject_id '_' teleporter '_Epochs_Entry.mat'],'file')
    question = 'You already have a list of epochs. What would you like to do?';
    questTitle = 'Redo pre-processing?';
    button1 = 'Use existing files';
    button2 = 'Redo pre-processing';
    answer = questdlg(question, questTitle, button1, button2, button2);
    
    if strcmpi(button2,answer) == 1
        redoEpochs = 1;
    end
end

if redoEpochs == 1
    
    % load the unity text file
    fid  = fopen(unityDataPath);
    
    if fixFile == 1
        txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n');
    else
        txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
    end
    
    fclose(fid);
    
    systemTime = txtdata{3};
    target     = txtdata{4};
    spaceType  = txtdata{6};
    timeType   = txtdata{7};
    xPos       = txtdata{8};
    zPos       = txtdata{9};
    yRot       = txtdata{10};
    
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
    
    % extract systemTime, space, and time, for each event
    teleporterEntryInd = logical(teleporterEntryInd);
    epochs = systemTime(teleporterEntryInd);
    epochsSpace = spaceType(teleporterEntryInd);
    epochsTime = timeType(teleporterEntryInd);
    
    % Convert epochs from ticks to EEG bins
    
    teleporterInds = find(teleporterEntryInd);
    epochsEDF1 = [];
    epochsEDF2 = [];
    missingEpochs = [];
    missingEpochInds = [];
    
    %  load tick/EEG conversions for each EDF
    load(['Mat Files/' save_stems{1} '_pulse_timing.mat']);
    indEEG1 = indEEG;
    unityTicks1 = unityTicks;
    
    load(['Mat Files/' save_stems{2} '_pulse_timing.mat']);
    indEEG2 = indEEG;
    unityTicks2 = unityTicks;
    
    % Loop through trials and find the onset time by fitting a line to the nearby pulses
    fitRange = 5;
    showFit = 0;
    for thisTrial = 1:length(teleporterInds)
        
        onsetTick = systemTime(teleporterInds(thisTrial));
        
        if (onsetTick - unityTicks1(end) < 0) % if trial in EDF1
            
            % find the pulse time (in ticks) closest to our tick
            tickDiff  = abs(onsetTick - unityTicks1);
            minDiffInd = find(tickDiff == min(tickDiff));
            
            % get a range of values before and after
            fitInds = [minDiffInd-fitRange:1:minDiffInd+fitRange];
            
            % get fit of line using the inds selected above
            fit_P = polyfit(unityTicks1(fitInds),indEEG1(fitInds),1);
            fit_y = polyval(fit_P,unityTicks1(fitInds));
            onsetBin = round(onsetTick*fit_P(1) + fit_P(2));
            
            if showFit == 1
                h = figure;
                plot(unityTicks1(fitInds),indEEG1(fitInds),'k*')
                hold on;
                plot(unityTicks1(fitInds),fit_y,'-')
                scatter(onsetTick,onsetBin,'ro')
                
                % make sure the fit looks good or else quit
                answer = questdlg('Does the fit look good?');
                if(strcmpi(answer,'Yes') ~= 1)
                    break
                end
                
                close(h);
            end
            
            epochsEDF1(end+1) = onsetBin;
            
        elseif (onsetTick - unityTicks2(1) < 0) % trial in lost EEG between EDF files
            missingEpochs(end + 1) = onsetTick;
            missingEpochInds(end + 1) = thisTrial;
            continue
        else % trial in EDF2
            
            % find the pulse time (in ticks) closest to our tick
            tickDiff  = abs(onsetTick - unityTicks2);
            minDiffInd = find(tickDiff == min(tickDiff));
            
            % get a range of values before and after
            fitInds = [minDiffInd-fitRange:1:minDiffInd+fitRange];
            
            % get fit of line using the inds selected above
            fit_P = polyfit(unityTicks2(fitInds),indEEG2(fitInds),1);
            fit_y = polyval(fit_P,unityTicks2(fitInds));
            onsetBin = round(onsetTick*fit_P(1) + fit_P(2));
            
            if showFit == 1
                h = figure;
                plot(unityTicks2(fitInds),indEEG2(fitInds),'k*')
                hold on;
                plot(unityTicks2(fitInds),fit_y,'-')
                scatter(onsetTick,onsetBin,'ro')
                
                % make sure the fit looks good or else quit
                answer = questdlg('Does the fit look good?');
                if(strcmpi(answer,'Yes') ~= 1)
                    break
                end
                
                close(h);
            end
            
            epochsEDF2(end+1) = onsetBin;
            
        end
    end
    
    
    % For each epoch, set the trial type
    
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
    epochs_saveFile = [subject_dir 'Mat Files/' subject_id '_' teleporter '_Epochs_Entry.mat'];
    save(epochs_saveFile,'epochsEDF1','epochsEDF2','missingEpochs','missingEpochInds','eSpace','eTime','eType');
end

%% Epoch the EEG data
load([subject_dir 'Mat Files/' subject_id '_' teleporter '_Epochs_Entry.mat']);
for thisDepth = 1:length(depthNames)
    
    % load the two EEG datasets
    clear EEG;
    EEG(1) = pop_loadset(['PreProcessing Intermediates/' save_stems{1} '_unepoched_' depthNames{thisDepth} '_marked.set']);
    EEG(2) = pop_loadset(['PreProcessing Intermediates/' save_stems{2} '_unepoched_' depthNames{thisDepth} '_marked.set']);
    
    % Insert events into EEG
    thisEpoch = 1;
    numSpikes = size(EEG(1).event,2);
    for n = numSpikes+1:numSpikes+length(epochsEDF1)
        
        EEG(1).event(n).latency = epochsEDF1(thisEpoch);
        
        if eSpace(thisEpoch) == 1 && eTime(thisEpoch) == 1
            EEG(1).event(n).type = '11';
        elseif eSpace(thisEpoch) == 1 && eTime(thisEpoch) == 2
            EEG(1).event(n).type = '12';
        elseif eSpace(thisEpoch) == 2 && eTime(thisEpoch) == 1
            EEG(1).event(n).type = '21';
        elseif eSpace(thisEpoch) == 2 && eTime(thisEpoch) == 2
            EEG(1).event(n).type = '22';
        end
        
        thisEpoch = thisEpoch + 1;
        
    end
    
    % Because of missing trials in the EEG recording break, identify which
    % epoch is the first of the EDF2 file
    eCount = length(eType) - length(epochsEDF2) + 1;
    
    % Insert events into EEG2
    thisEpoch = 1;
    numSpikes = size(EEG(2).event,2);
    for n = numSpikes+1:numSpikes+length(epochsEDF2)
        
        EEG(2).event(n).latency = epochsEDF2(thisEpoch);
        
        if eSpace(eCount) == 1 && eTime(eCount) == 1
            EEG(2).event(n).type = '11';
        elseif eSpace(eCount) == 1 && eTime(eCount) == 2
            EEG(2).event(n).type = '12';
        elseif eSpace(eCount) == 2 && eTime(eCount) == 1
            EEG(2).event(n).type = '21';
        elseif eSpace(eCount) == 2 && eTime(eCount) == 2
            EEG(2).event(n).type = '22';
        end
        
        eCount = eCount + 1;
        thisEpoch = thisEpoch + 1;
    end
    
    
    %% Save cleaned unepoched data
    for thisEDF = 1:length(save_stems)
        
        % Save unepoched data with spikes removed
        indexes = marks_label2index(EEG(thisEDF).marks.time_info, labelsToRemove_noSpikes, 'indexes','exact','on');
        
        if isempty(indexes)
            warning('No spikes found. Saving the data as it is.')
            noSpikeEEG(thisEDF) = pop_saveset(EEG(thisEDF), ['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisDepth} '_noSpikes.set']);
        else
            
            [noSpikeEEG(thisEDF),~] = pop_marks_select_data(EEG(thisEDF),'time marks',indexes);
            noSpikeEEG(thisEDF) = pop_saveset(noSpikeEEG(thisEDF), ['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisDepth} '_noSpikes.set']);
        end
        
        % if the spikes overlap with a real event, it will be removed, so check
        % for that now
        origEventLatency  = cell2mat({EEG(thisEDF).event.latency});
        [missingSpikeLatency, missingSpikeInd] = intersect(origEventLatency, indexes);
        if thisEDF == 1
            goodNoSpikeEpochs1 = [1:1:length(origEventLatency)]';
            goodNoSpikeEpochs1(missingSpikeInd) = [];
        else
            goodNoSpikeEpochs2 = [1:1:length(origEventLatency)]';
            goodNoSpikeEpochs2(missingSpikeInd) = [];
        end
        
        % Save version with spikes and waves removed
        indexes = marks_label2index(EEG(thisEDF).marks.time_info, labelsToRemove_noWaves, 'indexes', 'exact', 'on');
        
        if isempty(indexes)
            warning('No sharp waves found. Saving the data as it is.')
            noWaveEEG(thisEDF) = pop_saveset(noSpikeEEG(thisEDF), ['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisDepth} '_noSpikes_noWaves.set']);
        else
            [noWaveEEG(thisEDF),~] = pop_marks_select_data(EEG(thisEDF), 'time marks', indexes);
            noWaveEEG(thisEDF) = pop_saveset(noWaveEEG(thisEDF), ['PreProcessing Intermediates/' save_stems{thisEDF} '_unepoched_' depthNames{thisDepth} '_noSpikes_noWaves.set']);
        end
        
        % if the spikes or waves overlap with a real event, it will be removed,
        % so check for that now
        [missingWaveLatency, missingWaveInd] = intersect(origEventLatency, indexes);
        if thisEDF == 1
            goodNoWaveEpochs1 = [1:1:length(origEventLatency)]';
            goodNoWaveEpochs1(missingWaveInd) = [];
        else
            goodNoWaveEpochs2 = [1:1:length(origEventLatency)]';
            goodNoWaveEpochs2(missingWaveInd) = [];
        end
        
    end
    
    %% Save cleaned epoched data with spikes removed
    [newEEG(1), goodEpochs1] = pop_epoch(noSpikeEEG(1),{'11' '12' '21' '22'},[eStart eEnd]);
    [newEEG(2), goodEpochs2] = pop_epoch(noSpikeEEG(2),{'11' '12' '21' '22'},[eStart eEnd]);
    
    % Merge the two EEG files
    newEEG(3) = pop_mergeset(newEEG(1), newEEG(2));
    
    % Make a combined list of the epochs we kept
    epochOffset = length(eType) - length(epochsEDF2); % update goodEpochs2 indices so they match the true epoch number
    goodEpochs1 = goodNoSpikeEpochs1(goodEpochs1);
    goodEpochs2 = goodNoSpikeEpochs2(goodEpochs2);
    goodEpochs2 = goodEpochs2 + epochOffset;
    goodEpochs = cat(1, goodEpochs1, goodEpochs2);
    
    % Save epoched data
    pop_saveset(newEEG(3), ['Epoched Data/' subject_id '_' teleporter '_epoched_' depthNames{thisDepth} '_noSpikes.set']);
    
    % Save list of good epochs for this electrode
    save(['Mat Files/' subject_id '_' teleporter '_' depthNames{thisDepth} '_noSpikes_goodEpochs.mat'],'goodEpochs');
    
    %% Save cleaned epoched data with spikes AND waves removed
    [newEEG(1), goodEpochs1] = pop_epoch(noWaveEEG(1),{'11' '12' '21' '22'},[eStart eEnd]);
    [newEEG(2), goodEpochs2] = pop_epoch(noWaveEEG(2),{'11' '12' '21' '22'},[eStart eEnd]);
    
    % Merge the two EEG files
    newEEG(3) = pop_mergeset(newEEG(1), newEEG(2));
    
    % Make a combined list of the epochs we kept
    epochOffset = length(eType) - length(epochsEDF2); % update goodEpochs2 indices so they match the true epoch number
    goodEpochs1 = goodNoWaveEpochs1(goodEpochs1);
    goodEpochs2 = goodNoWaveEpochs2(goodEpochs2);
    goodEpochs2 = goodEpochs2 + epochOffset;
    goodEpochs = cat(1, goodEpochs1, goodEpochs2);
    
    % Save epoched data
    pop_saveset(newEEG(3), ['Epoched Data/' subject_id '_' teleporter '_epoched_' depthNames{thisDepth} '_noSpikes_noWaves.set']);
    
    % Save list of good epochs for this electrode
    save(['Mat Files/' subject_id '_' teleporter '_' depthNames{thisDepth} '_noSpikes_noWaves_goodEpochs.mat'],'goodEpochs');
    save(['Mat Files/' subject_id '_' teleporter '_' depthNames{thisDepth} '_noSpikes_noWaves_unepoched_goodEpochs.mat'], 'goodNoWaveEpochs1', 'goodNoWaveEpochs2');
    
end