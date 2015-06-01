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
        load([subjects{thisSubject} '_' versions{thisVersion} '_pepisode_by_condition.mat']);
        
        % average across conditions
        meanData = squeeze(mean(PepisodeByFreq,2));
        
        
        % plot the data for each electrode
        for thisChan = 1:size(meanData,2)
            
            thisData = squeeze(meanData(:,thisChan,:));
            h = figure;
            hold on;
            
            for thisCondition = 1:size(thisData,1)
                semilogx(F, thisData(thisCondition,:), 'LineWidth', 2, 'color', conditionColors{thisCondition});
                set(gca,...
                    'xtick', round(logspace(log10(F(1)),log10(F(end)),10)*100)/100,...
                    'xscale', 'log',...
                    'xlim', [F(1) F(end)],...
                    'ylim', [0 1]);
                
            end % thisCondition
            
            % plot frequency band boundaries
            ylims = get(gca,'ylim');
            plot([4 4],ylims,'--k','LineWidth',1);
            plot([8 8],ylims,'--k','LineWidth',1);
            plot([12 12],ylims,'--k','LineWidth',1);
            plot([30 30],ylims,'--k','LineWidth',1);
            hold off;
            
            legend(conditionNames)
            title(['Pepisode ' subjects{thisSubject} ' ' versions{thisVersion} ' ' chanList{thisChan}]);
            
            % save the figure
            save_dir = [exp_dir subjects{thisSubject} '/Figures/Pepisode_PSDs/'];
            if ~exist(save_dir, 'dir')
                system(['mkdir ' save_dir]);
            end
            
            fileName = ['Pepisode ' subjects{thisSubject} ' ' versions{thisVersion} ' ' chanList{thisChan} ' By TimePoint.png'];
            saveas(h, [save_dir fileName]);    
            
            close(h);
            
        end
        
    end % thisVersion
    
end % thisSubject


