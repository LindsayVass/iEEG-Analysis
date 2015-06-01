% Test whether pepisode varies within a frequency band based on timepoint,
% spatial condition, or temporal condition. This analysis will use the
% Kruskal-Wallis test (nonparametric one-way ANOVA) to
% determine whether pepisode significantly varies across conditions in a
% given frequency band for each electrode.
%
% Lindsay Vass 5 May 2015

clear all; close all; clc;

addpath('/Users/Lindsay/Documents/MATLAB/functions/');

subjects = {'UCDMC14' 'UCDMC15'};
versions = {'TeleporterA' 'TeleporterB'};

exp_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';

% where to save histogram
save_file = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/Hist_pepisode_by_time_space_timepoint.png';

% frequencies 
F = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011
deltaBand = find(F >= 1 & F <= 4);
thetaBand = find(F > 4 & F <= 8);
alphaBand = find(F > 8 & F <= 12);
betaBand  = find(F > 12 & F <= 30);
gammaBand = find(F > 30);
freqNames = {'deltaBand','thetaBand','alphaBand','betaBand','gammaBand'};

conditionNames = {'NSNT','NSFT','FSNT','FSFT'};


%% load in the data for each subject/version

pTime  = [];
pSpace = [];
pTimePoint = [];

for thisSubject = 1:length(subjects)
    
    for thisVersion = 1:length(versions)
        
        % navigate to subject's directory
        cd([exp_dir subjects{thisSubject} '/Mat Files/']);
        
        % load the mat file
        load([subjects{thisSubject} '_' versions{thisVersion} '_pepisode_all_epochs.mat']);
        
        for thisType = 1:size(allEpochData,2)
            
            dataCellHolder = allEpochData{thisType};
            
            allTimePointsData = [];
            for thisTimePoint = 1:size(dataCellHolder,1)
                
                tempData = dataCellHolder{thisTimePoint};
                allTimePointsData = cat(4, allTimePointsData, tempData);
                
            end % thisTimePoint
            
            % take the mean across time points
            thisTypeData = squeeze(nanmean(allTimePointsData,4));
            
            % Put data into a cell array that holds the epoch-wise data for
            % each condition
            conditionAnalysisCell(thisType) = {thisTypeData};
            
        end % thisType
        
        
        % Pull out the values for each frequency band / electrode and
        % compare across conditions
        pHolderSpace = nan(size(thisTypeData,1), length(freqNames));
        pHolderTime = pHolderSpace;
        pHolderTimePoint = pHolderSpace;
        for thisChan = 1:size(thisTypeData,1)
            
            for thisFreqBand = 1:length(freqNames)
                
                anovaData = [];
                anovaLabels = [];
                freqInd = eval(freqNames{thisFreqBand});
                
                for thisType = 1:size(allEpochData,2)
                    
                   tempData = conditionAnalysisCell{thisType};
                   tempData = squeeze(tempData(thisChan,:,freqInd));
                   
                   % take the mean across frequencies within this band
                   meanData = nanmean(tempData, 2);
                   
                   % add it to the summary array
                   anovaData = cat(1, anovaData, meanData);
                   anovaLabels = cat(1, anovaLabels, repmat(thisType,[size(meanData)]));
                    
                end % thisType
                
                % get indices for conditions of interest
                % TIME = NT (1,3) vs FT (2,4)
                % SPACE = NS (1,2) vs FS (3,4)
                NTind = find(anovaLabels == 1 | anovaLabels == 3);
                FTind = find(anovaLabels == 2 | anovaLabels == 4);
                NSind = find(anovaLabels == 1 | anovaLabels == 2);
                FSind = find(anovaLabels == 3 | anovaLabels == 4);
                
                timeData   = cat(1, anovaData(NTind), anovaData(FTind));
                timeLabels = cat(1, anovaLabels(NTind), anovaData(FTind));
                
                spaceData   = cat(1, anovaData(NSind), anovaData(FSind));
                spaceLabels = cat(1, anovaLabels(NSind), anovaLabels(FSind));
                
                % concatenate data for timepoint analysis
                timePointData = [];
                timePointLabels = [];
                for thisTimePoint = 1:size(allTimePointsData,4)
                    tempData = squeeze(allTimePointsData(thisChan, :, freqInd, thisTimePoint));
                    tempData = squeeze(mean(tempData, 2));
                    
                    timePointData = cat(1, timePointData, tempData);
                    timePointLabels = cat(1, timePointLabels, repmat(thisTimePoint, [size(allTimePointsData,2), 1]));
                end % thisTimePoint
                
                % perform kruskal-wallis test
                pHolderTime(thisChan, thisFreqBand)  = kruskalwallis(timeData, timeLabels, 'off');
                pHolderSpace(thisChan, thisFreqBand) = kruskalwallis(spaceData, spaceLabels, 'off');
                pHolderTimePoint(thisChan, thisFreqBand) = kruskalwallis(timePointData, timePointLabels, 'off');
                
            end % thisFreqBand
            
        end % thisChan
        
        % add to the big summary array
        pTime  = cat(1, pTime, pHolderTime);
        pSpace = cat(1, pSpace, pHolderSpace);
        pTimePoint = cat(1, pTimePoint, pHolderTimePoint);
        
    end % thisVersion
    
end % thisSubject

% find indices of p values < 0.05
pTimeSig = pTime < 0.05;
pSpaceSig = pSpace < 0.05;
pTimePointSig = pTimePoint < 0.05;

pTimeByFreqBand = mean(pTimeSig, 1);
pSpaceByFreqBand = mean(pSpaceSig, 1);
pTimePointByFreqBand = mean(pTimePointSig, 1);

%% plot results

% combine all into one matrix
allPData = cat(2, pTimeByFreqBand', pSpaceByFreqBand', pTimePointByFreqBand');

figure(1);
bar(allPData);
hold on;
plot(xlim, [0.05 0.05], '--k');
hold off;
title('PEPISODE: Proportion of electrodes modulated by space, time, and timepoint')
legend({'Time','Space','TimePoint'});
set(gca, ...
    'ylim', [0 0.5], ...
    'XTickLabel', freqNames)

saveas(1, save_file);
%% 

% binomial test on frequency of significant electrodes
myBinomTest(8, 32, 0.05, 'Greater')

