function Teleporter_Epoched_Pepisode(...
    subjectID, ...
    subjectDir, ...
    teleporter, ...
    chanList, ...
    cleanedUnepochedPrefix, ...
    cleanedUnepochedSuffix, ...
    cleanedEpochedPrefix,...
    cleanedEpochedSuffix,...
    saveFile, ...
    timePointNames, ...
    timesNT, ...
    timesFT, ...
    frequencies)
% function Teleporter_Epoched_Pepisode(...
%     subjectID, ...
%     subjectDir, ...
%     teleporter, ...
%     chanList, ...
%     cleanedUnepochedPrefix, ...
%     cleanedUnepochedSuffix, ...
%     cleanedEpochedPrefix,...
%     cleanedEpochedSuffix,...
%     saveFile, ...
%     timePointNames, ...
%     timesNT, ...
%     timesFT, ...
%     frequencies)
%
% EXAMPLE INPUTS:
%
% subjectID  = 'UCDMC13';
% subjectDir = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/'];
% teleporter = 'TeleporterA';
% chanList  = {'LAD1' 'LHD1'};
%
% % Specify file naming conventions for data. There are
% % separate EEG files for each depth electrode, so we will specify every
% % part of the path except the depth electrode name. Then, we will later
% % combine these as [prefix depthName suffix]. We will use a cell array to
% % allow for multiple prefixes or suffixes (e.g., different prefixes for
% % EDF1 and EDF2)
% cleanedUnepochedPrefix = {[subjectDir 'PreProcessing Intermediates/' subjectID '_' teleporter '_unepoched_']};
% cleanedUnepochedSuffix = {'_noSpikes_noWaves.set'};
% cleanedEpochedPrefix   = {[subjectDir 'Epoched Data/' subjectID '_' teleporter '_epoched_']};
% cleanedEpochedSuffix   = {'_noSpikes_noWaves.set'};
%
% % Specify path to save the cell array of pepisode values to. Don't
% include a file type because we'll save both csv and mat
% saveFile = [subjectDir 'Mat Files/Pepisode/Summary/' subjectID '_' teleporter '_epoched_pepisode_summary_noSpikes_noWaves_3sBuffer'];
%
% % time periods of interest in ms relative to teleporter entry
% timePointNames = {'Pre3' 'Pre2' 'Pre1' 'Tele' 'Post1' 'Post2' 'Post3'};
%
% timesNT = [-3000 -2001; -2000 -1001; -1000 0; 1 1830; 1831 2830; 2831 3830; 3831 4830];
% timesFT = [-3000 -2001; -2000 -1001; -1000 0; 1 2830; 2831 3830; 3831 4830; 4831 5830];
%
% % frequencies to use
% frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% set paths and filenames
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/iEEG/functions/'));

% Make the output directory if it doesn't already exist
if ~exist([subjectDir 'Mat Files/Pepisode/'], 'dir')
    system(['mkdir ' subjectDir 'Mat\ Files/Pepisode/']);
end

if ~exist([subjectDir 'Mat Files/Pepisode/Summary/'], 'dir')
    system(['mkdir ' subjectDir 'Mat\ Files/Pepisode/Summary/']);
end

if ~exist([subjectDir 'Mat Files/Pepisode/Power Distributions/'], 'dir')
    system(['mkdir ' subjectDir 'Mat\ Files/Pepisode/Power\ Distributions/']);
end

% Initialize the cell array to hold all of our pepisode values
pepisodeSummary      = cell(1,10);
pepisodeSummary(1,:) = {'SubjectID','Teleporter','Electrode','TrialNumber','TrialSpaceType','TrialTimeType','TrialType','TimePoint','Frequency','Pepisode'};
thisRow = 2;

% Get depth names from chanList
depthNames = unique(cellfun(@(s) s(1:3), chanList, 'UniformOutput', false));

%% Perform analysis
for thisDepth = 1:length(depthNames)
    
    % make a cell array of the EEG dataset paths
    eegList = cell(length(cleanedUnepochedPrefix), 1);
    for thisEEG = 1:length(cleanedUnepochedPrefix)
        eegPath = [cleanedUnepochedPrefix{thisEEG} depthNames{thisDepth} cleanedUnepochedSuffix{thisEEG}];
        eegList(thisEEG) = {eegPath};
    end
    
    % Find the channels on this depth electrode
    chanNames = findChansOnThisElectrode(chanList, depthNames{thisDepth});
    
    for thisChan = 1:length(chanNames)
        
        % Calculate power distribution if it doesn't already exist
        powerDistSaveFile = [subjectDir 'Mat Files/Pepisode/Power Distributions/' subjectID '_' teleporter '_' chanNames{thisChan} '_power_distribution.mat'];
        if ~exist(powerDistSaveFile, 'file')
            calcPowerDistributionLKV(eegList, chanNames{thisChan}, frequencies, 6, powerDistSaveFile);
        else
            load(powerDistSaveFile);
        end
        
        % Load the epoched EEG data
        epochedEEGPath = [cleanedEpochedPrefix{1} depthNames{thisDepth} cleanedEpochedSuffix{1}];
        EEG = pop_loadset(epochedEEGPath);
        
        % Find the index of this channel in the EEG.data
        chanLabels = {EEG.chanlocs.labels};
        chanInd = find(strcmpi(chanNames{thisChan}, chanLabels));
        
        % Create tables of epoch onsets and offsets
        trialTimesTable = makeTrialTimesTable(EEG, timesNT, timesFT, '1', '2');
        
        % Loop through trials
        for thisTrial = 1:size(trialTimesTable, 1)
            
            trialOnsets  = trialTimesTable{thisTrial, 2};
            trialOffsets = trialTimesTable{thisTrial, 3};
            
            % Loop through epochs within trial
            for thisTimePoint = 1:length(trialOnsets)
                
                timeInds = find(EEG.times >= trialOnsets(thisTimePoint) & EEG.times <= trialOffsets(thisTimePoint));
                eegData = EEG.data(chanInd, timeInds, thisTrial);
                
                [~, percentTimePepisode] = calcEpochedPepisodeLKV(powerDistribution, frequencies, eegData, EEG.srate, 95, 3);
                
                % Add to output array
                if strcmpi('1', EEG.event(thisTrial).type(1)) == 1
                    spaceType = 'NS';
                else
                    spaceType = 'FS';
                end
                if strcmpi('1', EEG.event(thisTrial).type(2)) == 1
                    timeType = 'NT';
                else
                    timeType = 'FT';
                end
                
                for thisFreq = 1:size(percentTimePepisode, 1)
                    pepisodeSummary(thisRow, :) = {subjectID, ...
                        teleporter, ...
                        chanNames{thisChan}, ...
                        thisTrial, ...
                        spaceType, ...
                        timeType, ...
                        [spaceType timeType], ...
                        timePointNames{thisTimePoint}, ...
                        frequencies(thisFreq), ...
                        percentTimePepisode(thisFreq)
                        };
                    thisRow = thisRow + 1;
                    
                end % thisFreq
                                               
                
            end % thisTimePoint
            
        end % thisTrial
        
    end % thisChan
    
end % thisDepth


% Write out the summary cell array to file
fprintf('\n\nWriting output to file...\n');
dlmcell([saveFile '.csv'], pepisodeSummary, 'delimiter', ',');
save([saveFile '.mat'], 'powerSummary');
