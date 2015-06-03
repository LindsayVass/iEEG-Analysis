% This script will take as input the file that was output by
% 'UCDMCXX_TeleporterX_Extract_Power' script and re-organize it to be in a
% tidy tall format for use in R. It will save the output as a csv.
%
% Lindsay Vass 6 May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

addpath('/Users/Lindsay/Documents/MATLAB/functions/');

subjects = {'UCDMC14' 'UCDMC15'};
versions = {'TeleporterA' 'TeleporterB'};

exp_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';

rOutFile = ['/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Power/Tall_Power_by_condition_' date '.csv'];

%% load in the data for each subject/version

allData = cell(length(subjects), length(versions));
for thisSubject = 1:length(subjects)
    
    for thisVersion = 1:length(versions)
        
        % navigate to subject's directory
        cd([exp_dir subjects{thisSubject} '/Mat Files/']);
        
        % load the mat file
        load([subjects{thisSubject} '_' versions{thisVersion} '_power.mat']);
        
        % put it into the holder array
        allData(thisSubject,thisVersion) = {allEpochData};
        
    end % thisVersion
    
end % thisSubject

%% organize data for R output

deltaBand = find(frequencies >= 1 & frequencies <= 4);
thetaBand = find(frequencies > 4 & frequencies <= 8);
alphaBand = find(frequencies > 8 & frequencies <= 12);
betaBand  = find(frequencies > 12 & frequencies <= 30);
gammaBand = find(frequencies > 30);
freqBandNames = {'deltaBand','thetaBand','alphaBand','betaBand','gammaBand'};
timePointNames = {'Pre' 'Tele' 'Post'};
conditionNames = {'NSNT','NSFT','FSNT','FSFT'};


rTable = {};
rTable(1,:) = {'Electrode','TrialType','FreqBand','TimePoint','Power'};
thisRow = 2;
varNames = {};
chanNumber = 1;
for thisDataset = 1:size(allData,1) * size(allData,2)
    
    thisData = allData{thisDataset};
    
    for thisType = 1:length(trialTypeList)
        
        thisTypeData = thisData{thisType};
        
        
        
        for thisTimePoint = 1:size(thisTypeData, 1)
            
            thisTimePointData = thisTypeData{thisTimePoint};
            
            % take the mean across epochs
            meanData = squeeze(mean(thisTimePointData, 2));
            
            for thisChan = 1:size(meanData,1)
                
                for thisFreq = 1:length(freqBandNames)
                    
                    % get the mean data across frequencies within this band
                    thisFreqData = squeeze(mean(meanData(thisChan, eval(freqBandNames{thisFreq})), 2));
                    
                    rTable(thisRow, :) = {num2str(chanNumber + thisChan - 1), trialTypeList{thisType}, freqBandNames{thisFreq}, timePointNames{thisTimePoint}, thisFreqData};
                    
                    thisRow = thisRow + 1;
                    
                end % thisFreq
                
            end % thisChan
            
        end % thisTimePoint
        
    end % thisType
    
    chanNumber = chanNumber + size(meanData,1);
    
end % thisDataset

% save the data
dlmcell(rOutFile, rTable, 'delimiter',',');
