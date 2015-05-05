% This script will perform a pepisode analysis, which consists of two
% parts. First, the script will calculate the power in each frequency
% across the entire experiment and use this to build a distribution of
% power values for each frequency. The value at 95% of the distribution is
% used as a power threshold in the second part of the analysis. For each
% epoch, the script will calculate the percent of time spent oscillating at
% a given frequency, where the oscillations must exceed both a power
% threshold (95% of the power distribution) and a duration threshold (# of
% cycles).
%
% In this script, we will evaluate three epochs:
% 1. Pre-teleportation (-1000 : 0 ms)
% 2. Teleportation (0 : 1830 ms for NT or 0 : 2830 ms for FT)
% 3. Post-teleportation (1830 : 1930 ms for NT or 2830 : 2930 ms for FT)
%
% Lindsay Vass 29 April 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% set parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/PepisodeCode/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/arne_code/'));
eeglab;

% set paths
subject_dir   = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/';
unepochedEEG1 = [subject_dir 'Raw Data/UCDMC15_TeleporterA_EDF1_unepoched.set']; % use the most pre-processed, but unepoched dataset
unepochedEEG2 = [subject_dir 'Raw Data/UCDMC15_TeleporterA_EDF2_unepoched.set']; % leave blank if only 1 EDF
epochsFile    = [subject_dir 'Mat Files/UCDMC15_TeleporterA_Epochs_Entry.mat']; % contains the onset and type of each epoch
unityFile     = [subject_dir 'Behavioral Data/TeleporterA/s3_FindStore_TeleporterA_FIXED.txt']; % Unity output during navigation to find stores
save_stem     = [subject_dir 'Mat Files/UCDMC15_TeleporterA_pepisode'];

% select pulse timing file
PTB_pulse_file = []; % time synchronization file for pulses from psychtoolbox
unity_EDF1_pulse_file = [subject_dir 'Raw Data/UCDMC15_TeleporterA_EDF1_pulse_timing.mat']; % ticks/bins for pulses from unity
unity_EDF2_pulse_file = [subject_dir 'Raw Data/UCDMC15_TeleporterA_EDF2_pulse_timing.mat']; % leave blank if only 1 EDF

% channel names to use
chanList = {'RAD4' 'RAD5' 'RAD6' 'RHD2' 'RHD3' 'RHD4' 'LAD3' 'LAD4' 'LHD1' 'LHD2' 'LHD3'};

% trial type names and corresponding eType
trialTypeList = {'NSNT' 'NSFT' 'FSNT' 'FSFT'};
trialeTypeList = [1:1:4];

% thresholds for pepisode 
durationThresh  = 3; % # of cycles
amplitudeThresh = 95; % percent of distribution

% if the distribution of power across the recording has already been
% calculated, set this to 1
skipCompute = 1;

% frequencies for a bandpass filter, to filter out the 60 Hz line noise
filtFreq = [59.9 60.9];

% time periods of interest
preTele = [-1000 0];
teleNT  = [1 1830];
teleFT  = [1 2830];
postNT  = [1831 2830];
postFT  = [2831 3830];

% frequencies to use
F = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% load EEG data
[EEG1] = pop_loadset(unepochedEEG1);

% find the indices of the channels of interest
chans1 = [];
for thisChan = 1:size(EEG1.data,1)
    thisChanName = {EEG1.chanlocs(thisChan).labels};
    goodChanInd = strcmpi(thisChanName,chanList);
    if sum(goodChanInd(:)) > 0
        chans1(end+1) = thisChan;
    end
end

if ~isempty(unepochedEEG2)
    [EEG2] = pop_loadset(unepochedEEG2);
    
    % find the indices of the channels of interest
    chans2 = [];
    for thisChan = 1:size(EEG2.data,1)
        thisChanName = {EEG2.chanlocs(thisChan).labels};
        goodChanInd = strcmpi(thisChanName,chanList);
        if sum(goodChanInd(:)) > 0
            chans2(end+1) = thisChan;
        end
    end
end


%% load unity file and identify start and end times of navigation
fid  = fopen(unityFile);
% txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n'); % use this version for unity output that we fixed in Matlab
txtdata =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1); % use this version for raw unity output
fclose(fid);
 
systemTime = txtdata{3};
startTicks = systemTime(1);
endTicks   = systemTime(end);

% make sure there's only one kind of pulse file
if ~isempty(PTB_pulse_file) && ~isempty(unity_EDF1_pulse_file)
    errorMsg = 'Error! You can only have one kind of pulse file. Get rid of either "PTB_pulse_file" or "unity_EDF1_pulse_file".';
    error(errorMsg);
end

% PTB pulses
if ~isempty(PTB_pulse_file)
    load(PTB_pulse_file);
    startBin = round(startTicks * time_sync_regression(1) + time_sync_regression(2));
    endBin   = round(endTicks * time_sync_regression(1) + time_sync_regression(2));
    
    % trim EEG data based on start and end times
    data1 = EEG1.data(chans1, startBin:endBin);
end

% unity pulses
bothTicks = cat(1,startTicks,endTicks);
if ~isempty(unity_EDF1_pulse_file)
    
    if ~isempty(unity_EDF2_pulse_file) % if 2 EDFs
        [~,EEG1bins,EEG2bins] = unityPulseTiming(bothTicks,unity_EDF1_pulse_file,unity_EDF2_pulse_file);
        
        if ~isempty(EEG1bins) && ~isempty(EEG2bins) % if start in EDF1 and end in EDF2
            startBin = EEG1bins(1);
            endBin   = EEG2bins(1);
            
            % trim EEG data
            data1 = EEG1.data(chans1, startBin:end);
            data2 = EEG2.data(chans2, 1:endBin);
            
        elseif ~isempty(EEG1bins) && isempty(EEG2bins) % if both in EDF1
            startBin = EEG1bins(1);
            endBin   = EEG1bins(2);
            
            % trim EEG data
            data1 = EEG1.data(chans1, startBin:endBin);
            
        elseif isempty(EEG1bins) && ~isempty(EEG2bins) % if both in EDF2
            startBin = EEG2bins(1);
            endBin   = EEG2bins(2);
            
            % trim EEG data
            data1 = EEG2.data(chans2, startBin:endBin);
        end
        
    else % if 1 EDF
        [~,EEG1bins,~] = unityPulseTiming(bothTicks,unity_EDF1_pulse_file);
        startBin = EEG1bins(1);
        endBin   = EEG1bins(2);
    end

end


%% Calculate pepisode for entire navigation period (not including free exploration)

%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute pepisode for the whole EEG file (this will return a
% vector with a 0 for every point in time without a sustained
% oscillation at that frequency, and a 1 for every point in time
% during which there is a sustained oscillation)

if ~skipCompute
    for thisChan = 1:size(data1,1)
        fprintf(['computing Pepisode for channel ' num2str(thisChan) ' of ' num2str(size(data1,1)) '\n']);
        eegdata = squeeze(data1(thisChan,:));
        outFile = strcat(save_stem,'_',EEG1.chanlocs(chans1(thisChan)).labels,'.mat');
        
        if exist('data2', 'var')
            eegdata2 = squeeze(data2(thisChan,:));
            calcPepisodeTwoEEG(eegdata,eegdata2,EEG1.srate,outFile,filtFreq,F,durationThresh,amplitudeThresh);
        else
            calcPepisode(eegdata,EEG1.srate,outFile,filtFreq,F,durationThresh,amplitudeThresh);
        end
    end
end

%% load epochs info
load(epochsFile);
if ~exist('epochsEDF1', 'var') && exist('epochs', 'var')
    epochsEDF1 = epochs';
    epochsEDF1   = round(epochsEDF1 * time_sync_regression(1) + time_sync_regression(2));
end

% get rid of rejected epochs
numEpochs1 = size(epochsEDF1,2);
eType1 = eType(1:numEpochs1);
epochsEDF1(rejectedEpochs1) = [];
eType1(rejectedEpochs1) = [];

% find onset times
onsets1 = cell(length(trialTypeList), 1);
for thisType = 1:length(trialTypeList)
    clear onsetList;
    onsetList = epochsEDF1(find(eType1 == trialeTypeList(thisType)));
    onsets1(thisType) = {onsetList};
end

if exist('epochsEDF2', 'var')
    
    % get rid of rejected epochs
    numEpochs2 = size(epochsEDF2,2);
    eType2 = eType(length(eType)-numEpochs2+1:end);
    epochsEDF2(rejectedEpochs2) = [];
    eType2(rejectedEpochs2) = [];
    
    % find onset times
    onsets2 = cell(length(trialTypeList), 1);
    for thisType = 1:length(trialTypeList)
        clear onsetList;
        onsetList = epochsEDF2(find(eType2 == trialeTypeList(thisType)));
        onsets2(thisType) = {onsetList};
    end
end

%%%%%%% step 2 %%%%%%%%%%%
% the second step in the Pepisode analysis is to obtain pepisode
% for every event of interest. The times are
% presented in terms of offset in samples in the eeg data file.

fprintf('\n\nretrieving Pepisode');

% parameters for retrieving Pepisode
intervalsNT = cat(1, preTele, teleNT, postNT);
durationsNT = round((intervalsNT(:,2) - intervalsNT(:,1))*EEG1.srate/1000);
intervalsFT = cat(1, preTele, teleFT, postFT);
durationsFT = round((intervalsFT(:,2) - intervalsFT(:,1))*EEG1.srate/1000);

% initialize summary arrays
PRE_1_all  = nan(length(trialTypeList), size(data1,1), size(F,2));
TELE_1_all = PRE_1_all;
POST_1_all = PRE_1_all;
PRE_2_all  = PRE_1_all;
TELE_2_all = PRE_1_all;
POST_2_all = PRE_1_all;

allEpochDataDimNames = {'Electrodes x Epochs x Frequencies'};

for thisType = 1:length(trialTypeList) % loop through trial types
    
    % if NT
    if ~isempty(strfind('NT', trialTypeList{thisType}))
        PREunionVecHolder1 = nan(size(data1,1), size(onsets1{thisType}, 2), size(F,2), durationsNT(1));
        TELEunionVecHolder1 = nan(size(data1,1), size(onsets1{thisType}, 2), size(F,2), durationsNT(2));
        POSTunionVecHolder1 = nan(size(data1,1), size(onsets1{thisType}, 2), size(F,2), durationsNT(3));
        
        if exist('onsets2', 'var')
            PREunionVecHolder2 = nan(size(data1,1), size(onsets2{thisType}, 2), size(F,2), durationsNT(1));
            TELEunionVecHolder2 = nan(size(data1,1), size(onsets2{thisType}, 2), size(F,2), durationsNT(2));
            POSTunionVecHolder2 = nan(size(data1,1), size(onsets2{thisType}, 2), size(F,2), durationsNT(3));
        end
        
        intervals = intervalsNT;
    else
        % if FT
        PREunionVecHolder1 = nan(size(data1,1), size(onsets1{thisType}, 2), size(F,2), durationsFT(1));
        TELEunionVecHolder1 = nan(size(data1,1), size(onsets1{thisType}, 2), size(F,2), durationsFT(2));
        POSTunionVecHolder1 = nan(size(data1,1), size(onsets1{thisType}, 2), size(F,2), durationsFT(3));
        
        if exist('onsets2', 'var')
            PREunionVecHolder2 = nan(size(data1,1), size(onsets2{thisType}, 2), size(F,2), durationsFT(1));
            TELEunionVecHolder2 = nan(size(data1,1), size(onsets2{thisType}, 2), size(F,2), durationsFT(2));
            POSTunionVecHolder2 = nan(size(data1,1), size(onsets2{thisType}, 2), size(F,2), durationsFT(3));
        end
        
        intervals = intervalsFT;
    end
    
    for thisChan = 1:size(data1,1); % loop through electrodes
        
        % file containing pepisode power info for each frequency
        outFile = strcat(save_stem,'_', chanList{thisChan},'.mat');
        if exist('onsets2', 'var')
            PREunionVecHolder1(thisChan,:,:,:)  = getpepisodeTwoEEG(outFile, 1, onsets1{thisType} - startBin, intervals(1,2) - intervals(1,1), intervals(1), F);
            TELEunionVecHolder1(thisChan,:,:,:) = getpepisodeTwoEEG(outFile, 1, onsets1{thisType} - startBin, intervals(2,2) - intervals(2,1), intervals(2), F);
            POSTunionVecHolder1(thisChan,:,:,:) = getpepisodeTwoEEG(outFile, 1, onsets1{thisType} - startBin, intervals(3,2) - intervals(3,1), intervals(3), F);
            
            PREunionVecHolder2(thisChan,:,:,:)  = getpepisodeTwoEEG(outFile, 2, onsets2{thisType}, intervals(1,2) - intervals(1,1), intervals(1), F);
            TELEunionVecHolder2(thisChan,:,:,:) = getpepisodeTwoEEG(outFile, 2, onsets2{thisType}, intervals(2,2) - intervals(2,1), intervals(2), F);
            POSTunionVecHolder2(thisChan,:,:,:) = getpepisodeTwoEEG(outFile, 2, onsets2{thisType}, intervals(3,2) - intervals(3,1), intervals(3), F);
            
        else
            PREunionVecHolder1(thisChan,:,:,:)  = getpepisode(outFile, onsets1{thisType} - startBin, intervals(1,2) - intervals(1,1), intervals(1), F);
            TELEunionVecHolder1(thisChan,:,:,:) = getpepisode(outFile, onsets1{thisType} - startBin, intervals(2,2) - intervals(2,1), intervals(2), F);
            POSTunionVecHolder1(thisChan,:,:,:) = getpepisode(outFile, onsets1{thisType} - startBin, intervals(3,2) - intervals(3,1), intervals(3), F);
        end
        
    end % thisChan
    
    % take the man across time so we can save out the data for each epoch
    mean_PRE_1_time = squeeze(nanmean(PREunionVecHolder1,4));
    mean_TELE_1_time = squeeze(nanmean(TELEunionVecHolder1,4));
    mean_POST_1_time = squeeze(nanmean(POSTunionVecHolder1,4));
    
    if exist('onsets2','var')
        mean_PRE_2_time = squeeze(nanmean(PREunionVecHolder2,4));
        mean_TELE_2_time = squeeze(nanmean(TELEunionVecHolder2,4));
        mean_POST_2_time = squeeze(nanmean(POSTunionVecHolder2,4));
        
        % if there aren't any trials in one set, use the other; otherwise
        % concatenate
        if length(size(mean_PRE_1_time)) == 2
            
            mean_PRE_1_time = mean_PRE_2_time;
            mean_TELE_1_time = mean_TELE_2_time;
            mean_POST_1_time = mean_POST_2_time;
            
        elseif length(size(mean_PRE_2_time)) == 2
        else
            mean_PRE_1_time = cat(2, mean_PRE_1_time, mean_PRE_2_time);
            mean_TELE_1_time = cat(2, mean_TELE_1_time, mean_TELE_2_time);
            mean_POST_1_time = cat(2, mean_POST_1_time, mean_POST_2_time);
        end
    end
    
    tempdata = {mean_PRE_1_time; mean_TELE_1_time; mean_POST_1_time};
    allEpochData(thisType) = {tempdata};
    
    % take the mean across epochs
    mean_PRE_1  = squeeze(nanmean(PREunionVecHolder1,2));
    mean_TELE_1 = squeeze(nanmean(TELEunionVecHolder1,2));
    mean_POST_1 = squeeze(nanmean(POSTunionVecHolder1,2));
    
    if exist('onsets2', 'var')
        mean_PRE_2  = squeeze(nanmean(PREunionVecHolder2,2));
        mean_TELE_2 = squeeze(nanmean(TELEunionVecHolder2,2));
        mean_POST_2 = squeeze(nanmean(POSTunionVecHolder2,2));
    end
    
    % take the mean across epochs
    mean_PRE_1  = squeeze(nanmean(mean_PRE_1,3));
    mean_TELE_1 = squeeze(nanmean(mean_TELE_1,3));
    mean_POST_1 = squeeze(nanmean(mean_POST_1,3));
    
    if exist('onsets2', 'var')
        mean_PRE_2  = squeeze(nanmean(mean_PRE_2,3));
        mean_TELE_2 = squeeze(nanmean(mean_TELE_2,3));
        mean_POST_2 = squeeze(nanmean(mean_POST_2,3));
    end
    
    % put into summary arrays
    PRE_1_all(thisType,:,:)  = mean_PRE_1;
    TELE_1_all(thisType,:,:) = mean_TELE_1;
    POST_1_all(thisType,:,:) = mean_POST_1;
    
    if exist('onsets2', 'var')
        PRE_2_all(thisType,:,:)  = mean_PRE_2;
        TELE_2_all(thisType,:,:) = mean_TELE_2;
        POST_2_all(thisType,:,:) = mean_POST_2;
    end
    
    
end % thisType

% save data for all epochs
save([save_stem '_all_epochs.mat'], 'allEpochData','allEpochDataDimNames');

% find nans
PRE_1_nan  = isnan(PRE_1_all);
TELE_1_nan = isnan(TELE_1_all);
POST_1_nan = isnan(POST_1_all);

if exist('onsets2' ,'var')
    PRE_2_nan  = isnan(PRE_2_all);
    TELE_2_nan = isnan(TELE_2_all);
    POST_2_nan = isnan(POST_2_all);
end

% if one condition is all NaN, this means there were no trials of that type
% in that EEG, so use the other EEG's data instead
if ~exist('onsets2', 'var')
    PRE_all = PRE_1_all;
    TELE_all = TELE_1_all;
    POST_all = POST_1_all;
else
    
    for thisType = 1:length(trialTypeList)
        
        % if EEG1 missing data
        if sum(sum(PRE_1_nan(thisType,:,:))) == size(PRE_1_nan,2) * size(PRE_1_nan,3)
            PRE_all(thisType,:,:)  = PRE_2_all(thisType,:,:);
            TELE_all(thisType,:,:) = TELE_2_all(thisType,:,:);
            POST_all(thisType,:,:) = POST_2_all(thisType,:,:);
            % if EEG2 missing data
        elseif sum(sum(PRE_2_nan(thisType,:,:))) == size(PRE_2_nan,2) * size(PRE_2_nan,3)
            PRE_all(thisType,:,:) = PRE_1_all(thisType,:,:);
            TELE_all(thisType,:,:) = TELE_1_all(thisType,:,:);
            POST_all(thisType,:,:) = POST_1_all(thisType,:,:);
        else % both have data
            PRE_all(thisType,:,:) = (PRE_1_all(thisType,:,:) + PRE_2_all(thisType,:,:)) / 2;
            TELE_all(thisType,:,:) = (TELE_1_all(thisType,:,:) + TELE_2_all(thisType,:,:)) / 2;
            POST_all(thisType,:,:) = (POST_1_all(thisType,:,:) + POST_2_all(thisType,:,:)) / 2;
            
        end
        
    end

end

% combine into one array
PepisodeByFreq = nan(3, size(PRE_all,1), size(PRE_all,2), size(PRE_all,3));
PepisodeByFreq(1,:,:,:) = PRE_all;
PepisodeByFreq(2,:,:,:) = TELE_all;
PepisodeByFreq(3,:,:,:) = POST_all;


%% get mean for each frequency band
deltaBand = find(F >= 1 & F <= 4);
thetaBand = find(F > 4 & F <= 8);
alphaBand = find(F > 8 & F <= 12);
betaBand  = find(F > 12 & F <= 30);
gammaBand = find(F > 30);
freqNames = {'deltaBand','thetaBand','alphaBand','betaBand','gammaBand'};

% initialize summary array; 3 (pre/tele/post) x Trial Type x Electrode x
% Frequency Band
PepisodeByFreqBand = nan(3, size(PRE_1_all,1), size(PRE_1_all,2), length(freqNames));

for thisType = 1:length(trialTypeList)
    for thisFreq = 1:length(freqNames)
        freqs = eval(freqNames{thisFreq});
        PRE_data = squeeze(PRE_all(thisType,:,freqs));
        TELE_data = squeeze(TELE_all(thisType,:,freqs));
        POST_data = squeeze(POST_all(thisType,:,freqs));
        
        % take the mean across frequencies
        PepisodeByFreqBand(1,thisType,:,thisFreq) = mean(PRE_data,2);
        PepisodeByFreqBand(2,thisType,:,thisFreq) = mean(TELE_data,2);
        PepisodeByFreqBand(3,thisType,:,thisFreq) = mean(POST_data,2);
        
    end % thisFreq
end  % thisType

% Reorder dimensions so channels x condition x freqband x time interval
% (pre/tele/post)
PepisodeByFreqBand = permute(PepisodeByFreqBand, [3 2 4 1]);

% save data
save([save_stem '_by_condition.mat'],'PepisodeByFreq','PepisodeByFreqBand','F','trialTypeList','chanList');

% %% plot mean across channels
meanPepisode = squeeze(mean(PepisodeByFreqBand,1));

h = figure;

ax = subplot(221);
bar(squeeze(meanPepisode(1, :, :)));
set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
title('NSNT')

ax = subplot(222);
bar(squeeze(meanPepisode(2, :, :)));
set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
title('NSFT')

ax = subplot(223);
bar(squeeze(meanPepisode(3, :, :)));
set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
title('FSNT')

ax = subplot(224);
bar(squeeze(meanPepisode(4, :, :)));
set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
title('FSFT')
% 
% 
% % make plots for each electrode
% for thisChan = 1:size(PepisodeByFreqBand,1)
%     
%     h = figure;
%     
%     ax = subplot(221);
%     bar(squeeze(PepisodeByFreqBand(thisChan, 1, :, :)));
%     set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
%     title('NSNT')
%     
%     ax = subplot(222);
%     bar(squeeze(PepisodeByFreqBand(thisChan, 2, :, :)));
%     set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
%     title('NSFT')
%     
%     ax = subplot(223);
%     bar(squeeze(PepisodeByFreqBand(thisChan, 3, :, :)));
%     set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
%     title('FSNT')
%     
%     ax = subplot(224);
%     bar(squeeze(PepisodeByFreqBand(thisChan, 4, :, :)));
%     set(gca,'XTickLabel',{'delta','theta','alpha','beta','gamma'});
%     title('FSFT')
% end
