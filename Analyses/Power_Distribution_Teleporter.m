% This script will calculate the power distribution for each electrode
% throughout the entire navigation period (excluding free exploration).
%
% Lindsay Vass 8 May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% set parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/PepisodeCode/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/arne_code/'));
eeglab;

% set paths
subject_dir   = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/';
unepochedEEG1 = [subject_dir 'Raw Data/UCDMC14_TeleporterA_unepoched.set']; % use the most pre-processed, but unepoched dataset
unepochedEEG2 = []; % leave blank if only 1 EDF
epochsFile    = [subject_dir 'Mat Files/UCDMC14_TeleporterA_Epochs_Entry.mat']; % contains the onset and type of each epoch
unityFile     = [subject_dir 'Behavioral Data/TeleporterA/s2_patientTeleporterData.txt']; % Unity output during navigation to find stores
saveStem     = [subject_dir 'Mat Files/UCDMC14_TeleporterA_'];

% select pulse timing file
PTB_pulse_file = [subject_dir 'Mat Files/UCDMC14_TeleporterA_time_sync.mat']; % time synchronization file for pulses from psychtoolbox
unity_EDF1_pulse_file = []; % ticks/bins for pulses from unity
unity_EDF2_pulse_file = []; % leave blank if only 1 EDF

% channel names to use
chanList = {'LAD1' 'LHD1' 'RAD1' 'RHD1' 'RHD2'};

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


%% Calculate power distribution for entire navigation period (not including free exploration)
for thisChan = 1:size(data1,1)
    fprintf(['computing power distribution for channel ' num2str(thisChan) ' of ' num2str(size(data1,1)) '\n']);
    eegdata = squeeze(data1(thisChan,:));
    
    if exist('data2', 'var')
        eegdata2 = squeeze(data2(thisChan,:));
        [powerDist, freqs] = calcPowerDistributionTwoEEG(eegdata,eegdata2,EEG1.srate,F,6);
    else
        [powerDist, freqs] = calcPowerDistribution(eegdata,EEG1.srate,F,6);
    end
    
    % save the distribution
    save([saveStem  'power_distribution_' chanList{thisChan} '.mat'], 'powerDist', 'freqs');
end

