close all;clear all; clc;

addpath(genpath('/Users/lindsay/Documents/MATLAB/eeglab13_3_2b'));
addpath(genpath('/Users/mcopara/Code/iEEG_code/code/downloaded_code/eeg_toolbox/eeg_toolbox_v1_3_2'));
addpath(genpath('/Users/mcopara/Code/iEEG_code/code/arne_code/matlab_scripts_from_ucla/'));

subjects = {'UCDMC13' 'UCDMC14'};

pre_NSNT = [];
pre_NSFT = [];
pre_FSNT = [];
pre_FSFT = [];

tele_NSNT = [];
tele_NSFT = [];
tele_FSNT = [];
tele_FSFT = [];

for thisSubject = 1:length(subjects)
    load(['/Users/lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/' subjects{thisSubject} '_pepisode.mat']);
    
    pre_NSNT = cat(2,pre_NSNT,mean_pre_NSNT);
    pre_NSFT = cat(2,pre_NSFT,mean_pre_NSFT);
    pre_FSNT = cat(2,pre_FSNT,mean_pre_FSNT);
    pre_FSFT = cat(2,pre_FSFT,mean_pre_FSFT);
    
    tele_NSNT = cat(2,tele_NSNT,mean_tele_NSNT);
    tele_NSFT = cat(2,tele_NSFT,mean_tele_NSFT);
    tele_FSNT = cat(2,tele_FSNT,mean_tele_FSNT);
    tele_FSFT = cat(2,tele_FSFT,mean_tele_FSFT);
    
end

% transpose so rows = electrodes and columns = frequency bands
pre_NSNT = pre_NSNT';
pre_NSFT = pre_NSFT';
pre_FSNT = pre_FSNT';
pre_FSFT = pre_FSFT';

tele_NSNT = tele_NSNT';
tele_NSFT = tele_NSFT';
tele_FSNT = tele_FSNT';
tele_FSFT = tele_FSFT';

T = table(pre_NSNT,tele_NSNT,pre_NSFT,tele_NSFT,pre_FSNT,tele_FSNT,pre_FSFT,tele_FSFT,'VariableNames',{'Pre_NSNT','Tele_NSNT','Pre_NSFT','Tele_NSFT','Pre_FSNT','Tele_FSNT','Pre_FSFT','Tele_FSFT'});

