% Make a line plot showing pepisode across frequencies for each condition
% for each electrode
%
% Lindsay Vass 5 May 2015

clear all; close all; clc;

subjects = {'UCDMC14' 'UCDMC15'};
versions = {'TeleporterA' 'TeleporterB'};

exp_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';

% Matlab colors
conditionColors = {'b','g','r'};
conditionNames = {'Pre','Tele','Post'};


%% load in the data for each subject/version

for thisSubject = 1:length(subjects)
    
    for thisVersion = 1:length(versions)
        
        % navigate to subject's directory
        cd([exp_dir subjects{thisSubject} '/Mat Files/']);
        
        % load the mat file
        load([subjects{thisSubject} '_' versions{thisVersion} '_power.mat']);
        
        % Loop through conditions and average across epochs
        preData  = nan(length(trialTypeList), length(chanList), length(frequencies));
        teleData = preData;
        postData = preData;
        for thisTrialType = 1:size(allEpochData, 2)
            thisTypeData = allEpochData{thisTrialType};
            
            tempPre = thisTypeData{1};
            tempPre = squeeze(mean(tempPre, 2));
            preData(thisTrialType,:,:) = tempPre;
            
            tempTele = thisTypeData{2};
            tempTele = squeeze(mean(tempTele, 2));
            teleData(thisTrialType,:,:) = tempTele;
            
            tempPost = thisTypeData{3};
            tempPost = squeeze(mean(tempPost, 2));
            postData(thisTrialType,:,:) = tempPost;
        end
        
        % concatenate the three time points
        allTimePoints = cat(4, preData, teleData, postData);
        allTimePoints = permute(allTimePoints, [4 1 2 3]);
        
        % average across conditions
        meanData = squeeze(mean(allTimePoints,2));
        
        
        % plot the data for each electrode
        for thisChan = 1:size(meanData,2)
            
            thisData = squeeze(meanData(:,thisChan,:));
            h = figure;
            hold on;
            
            for thisCondition = 1:size(thisData,1)
                loglog(frequencies, log(thisData(thisCondition,:)), 'LineWidth', 2, 'color', conditionColors{thisCondition});
                set(gca,...
                    'xtick', round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,...
                    'xscale', 'log',...
                    'xlim', [frequencies(1) frequencies(end)]);
                
            end % thisCondition
            
            % plot frequency band boundaries
            ylims = get(gca,'ylim');
            plot([4 4],ylims,'--k','LineWidth',1);
            plot([8 8],ylims,'--k','LineWidth',1);
            plot([12 12],ylims,'--k','LineWidth',1);
            plot([30 30],ylims,'--k','LineWidth',1);
            hold off;
            
            legend(conditionNames)
            title(['Power ' subjects{thisSubject} ' ' versions{thisVersion} ' ' chanList{thisChan}]);
            
            % save the figure
            save_dir = [exp_dir subjects{thisSubject} '/Figures/Power_PSDs/'];
            if ~exist(save_dir, 'dir')
                system(['mkdir ' save_dir]);
            end
            
            fileName = ['Power ' subjects{thisSubject} ' ' versions{thisVersion} ' ' chanList{thisChan} ' By TimePoint.png'];
            saveas(h, [save_dir fileName]);    
            
            close(h);
            
        end
        
    end % thisVersion
    
end % thisSubject


