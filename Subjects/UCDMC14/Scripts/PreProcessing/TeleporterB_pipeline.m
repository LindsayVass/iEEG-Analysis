% Pre-process UCDMC14 data using new pipeline

%% set up file naming
eeglab;

% experiment directory, contains a folder for each patient
expDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';

% patient ID, same as name of patient folder (expDir/subjID)
subjID  = 'UCDMC14';
subjDir = [expDir subjID '/'];
saveStem = 'UCDMC14_TeleporterB_';

%% load raw data 
edfPath = [subjDir 'Raw Data/UCDMC14_020715.edf'];
EEG     = edf2eeg(edfPath, subjID);

%% save raw data as EEG mat
preprocDir = [subjDir 'PreProcessing Intermediates/Pipeline/'];
if ~exist(preprocDir, 'dir')
    mkdir(preprocDir);
end

% save
savePath = [preprocDir saveStem 'raw.set'];
pop_saveset(EEG, savePath);

%% check for electrodes with gross signal artifacts
EEG = pop_loadset(savePath);
eeglab redraw;

% mark the marker channel and any obviously bad channels
EEG = markBadChannels(EEG);
keyboard;

% save updated version
savePath = [preprocDir saveStem 'C01.set'];
pop_saveset(EEG, savePath);

%% re-reference the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using all good channels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
EEGreref = rerefAllGoodChans(EEG);

% view original data and rereferenced data together
% original = blue
% re-ref = gray
% excluded channels = gray with gray/red triangle next to channel name
EEG = pop_vised(EEG, 'data2', 'EEGreref.data');

% save updated version
savePath = [preprocDir saveStem 'C01_Reref_All.set'];
pop_saveset(EEGreref, savePath);

%% perform artifact detection/removal

% prepare data
origDataPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/PreProcessing Intermediates/Pipeline/UCDMC14_TeleporterB_C01_Reref_All.set';
EEG = pop_loadset(origDataPath);
eeglab redraw;

outputDir = [preprocDir 'UCDMC14_TeleporterB_C01_Reref_All/'];
outputStem = 'UCDMC14_TeleporterB_C01_Reref_All_';

% Clean the data for each channel separately. This next function will first
% split your dataset into multiple data sets containing one channel each.
% These are contained in a folder called 'dirty_unepoched', which is
% created inside outputDir. It will then epoch each channel's data into
% short consecutive epochs (length = epochSecs). Finally, it will flag any
% epoch that contains extreme values, defined as a value that exceeds a
% certain number of standard deviations away from the mean (threshold =
% numSD). These marked and epoched files will be contained in a folder
% called 'clean_epoched', also found in outputDir. 
%
% If your data set does not divide evenly by epochSecs, it will pad the
% data with NaN and return the number of samples to later trim when the
% data is recombined (samplesToTrim).
epochSecs = 1;
numSD     = 7;
[splitFileList, markerPath, samplesToTrim] = splitAndCleanDataset(EEG, outputDir, outputStem, epochSecs, numSD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 2: Recombine channels on the same strip/grid/depth %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For this to work, channels on the same strip must have the same string
% (e.g., LAD1, LAD2, LAD3 all share 'LAD')           
mergeFileList = mergeDatasetsByStrip(splitFileList, samplesToTrim, outputDir, outputStem);

% view time points marked for rejection (do for each strip by changing the
% number in fileList{1} below
EEG = pop_loadset(mergeFileList{1});
EEG = pop_vised(EEG);

%% clip artifact events from merged data
outputDirNoSpace = strrep(outputDir, ' ', '\ ');
system(['mkdir ' outputDirNoSpace 'clean_merged/']);
for thisFile = 1:length(mergeFileList)
    EEG = pop_loadset(mergeFileList{thisFile});
    marks_ind = marks_label2index(EEG.marks.time_info, {'manual', 'rejthresh'}, 'indexes', 'exact', 'on');
    EEG = pop_marks_select_data(EEG, 'time marks', marks_ind);
    savePath = [outputDir 'clean_merged/' outputStem EEG.chanlocs(1).labels(1:3) '_clean.set'];
    pop_saveset(EEG, savePath);
end

%% add experiment events
epochsPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Mat Files/UCDMC14_TeleporterB_Epochs_Entry.mat';
load(epochsPath);
eStart = -3;
eEnd   = 6;

for thisFile = 1:length(mergeFileList)
    EEG = pop_loadset(mergeFileList{thisFile});
    for thisTrial = 1:length(epochsEDF)
        thisLabel = strcat(num2str(eSpace(thisTrial)), num2str(eTime(thisTrial)));
        EEG.event = addExperimentEvent(EEG.event, epochsEDF(thisTrial), thisLabel);
    end
    
    % save cleaned unepoched data
    ind = marks_label2index(EEG.marks.time_info, {'manual', 'rejthresh'}, 'indexes', 'exact', 'on');
    savePath = [mergeFileList{thisFile}(1:end-4) '_clean.set'];
    if isempty(ind)
        warning('No artifacts found. Saving the data as it is.')
        cleanEEG = EEG;
        pop_saveset(EEG, savePath);
    else
        [cleanEEG, ~] = pop_marks_select_data(EEG, 'time marks', ind);
        pop_saveset(cleanEEG, savePath);
    end
    
    % if artifacts overlap with experimental event, it will be removed
    % silently, so check for that now
    origLatency = cell2mat({EEG.event.latency});
    [missingLatency, missingInd] = intersect(origLatency, ind);
    noArtifactEpochs = [1:1:length(origLatency)]';
    noArtifactEpochs(missingInd) = [];
    
    % epoch data
    [epochEEG, goodEpochs] = pop_epoch(cleanEEG, {'11' '12' '21' '22'}, [eStart eEnd]);
    goodEpochs = noArtifactEpochs(goodEpochs);
    
    % save epoched data
    stripName = mergeFileList{thisFile}(end-6:end-4);
    epochSavePath = [subjDir 'Epoched Data/Pipeline/UCDMC14_TeleporterB_' stripName '_epoched.set'];
    pop_saveset(epochEEG, epochSavePath);
    
    matSavePath = [subjDir 'Mat Files/Pipeline/UCDMC14_TeleporterB_' stripName '_goodEpochs.mat'];
    save(matSavePath, 'goodEpochs');
    
end
