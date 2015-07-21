function PepisodeSustainedness(...
    subjectID, ...
    subjectDir, ...
    teleporter, ...
    chanList, ...
    cleanedUnepochedPrefix, ...
    cleanedUnepochedSuffix, ...
    cleanedEpochedPrefix,...
    cleanedEpochedSuffix,...
    analysisDir, ...
    saveFile, ...
    timePointNames, ...
    timesNT, ...
    timesFT, ...
    frequencies)
% function PepisodeSustainedness(...
%     subjectID, ...
%     subjectDir, ...
%     teleporter, ...
%     chanList, ...
%     cleanedUnepochedPrefix, ...
%     cleanedUnepochedSuffix, ...
%     cleanedEpochedPrefix,...
%     cleanedEpochedSuffix,...
%     analysisDir, ...
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
% % Specify path to save the csv of pepisode values to. Main
% directory is specified by analysisDir:
% analysisDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Epoched_Equal_Time_Bins/';
%
% Within that directory, the script will save files using the saveFile
% provided:
% saveFile = [subjectID '_' teleporter '_pepisode_sustainedness.csv'];
% It will save a .csv to analysisDir/csv/
%
% time periods of interest in ms relative to teleporter entry
% timePointNames = {'WholeTrial'};
% 
% timesNT = [-1830 3660];
% timesFT = [-2830 5660];
%
% % frequencies to use
% frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% set paths and filenames
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_4_4b'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/iEEG/functions/'));

% escape spaces in the directory name
analysisDirNoSpace = strrep(analysisDir, ' ', '\ ');

% Make the output directory if it doesn't already exist
if ~exist([analysisDir], 'dir')
    system(['mkdir ' analysisDirNoSpace]);
end

if ~exist([analysisDir 'csv/'], 'dir')
    system(['mkdir ' analysisDirNoSpace 'csv/']);
end


observationID = 1;

%% Get depth names from chanList
depthNames = unique(cellfun(@(s) s(1:3), chanList, 'UniformOutput', false));

%% Perform analysis
for thisDepth = 1:length(depthNames)
    
    fprintf(['Loading electrode ' depthNames{thisDepth} '\n\n']);
    
    % make a cell array of the EEG dataset paths
    eegList = cell(length(cleanedUnepochedPrefix), 1);
    for thisEEG = 1:length(cleanedUnepochedPrefix)
        eegPath = [cleanedUnepochedPrefix{thisEEG} depthNames{thisDepth} cleanedUnepochedSuffix{thisEEG}];
        eegList(thisEEG) = {eegPath};
    end
    
    % Find the channels on this depth electrode
    chanNames = findChansOnThisElectrode(chanList, depthNames{thisDepth});
    
    for thisChan = 1:length(chanNames)
        
        %% Initialize the output cell arrays
        
        % This array will hold the details about each observation
        observationCharacteristics      = cell(1,10);
        observationCharacteristics(1,:) = {'ObservationID','SubjectID','Teleporter','Electrode','TrialNumber','TrialSpaceType','TrialTimeType','TrialType','TimePoint','Frequency'};
        
        % our other output arrays will be dataNT and dataFT, but we won't know how
        % wide they are until we run the first subject's analysis, so we'll create
        % them then
        
        charRow = 2;
        ntRow = 2;
        ftRow = 2;
        
        % Calculate power distribution if it doesn't already exist
        powerDistSaveFile = [subjectDir 'Mat Files/Pepisode/Power Distributions/' subjectID '_' teleporter '_' chanNames{thisChan} '_power_distribution.mat'];
        if ~exist(powerDistSaveFile, 'file')
            powerDistribution = calcPowerDistributionLKV(eegList, chanNames{thisChan}, frequencies, 6, powerDistSaveFile);
        else
            load(powerDistSaveFile);
        end
        
        % Load the epoched EEG data
        epochedEEGPath = [cleanedEpochedPrefix{1} depthNames{thisDepth} cleanedEpochedSuffix{1}];
        EEG = pop_loadset(epochedEEGPath);
        
        fprintf(['\n\nChan ' num2str(thisChan) ': ']);
        
        % Find the index of this channel in the EEG.data
        chanLabels = {EEG.chanlocs.labels};
        chanInd = find(strcmpi(chanNames{thisChan}, chanLabels));
        
        % Create tables of epoch onsets and offsets
        trialTimesTable = makeTrialTimesTable(EEG, timesNT, timesFT, '1', '2'); 
        
        % Loop through trials
        for thisTrial = 1:size(trialTimesTable, 1)
            
            fprintf([' ' num2str(thisTrial)]);
            
            trialOnsets  = trialTimesTable{thisTrial, 2};
            trialOffsets = trialTimesTable{thisTrial, 3};
            
            % Loop through epochs within trial
            for thisTimePoint = 1:length(trialOnsets)
                
                timeInds = find(EEG.times >= trialOnsets(thisTimePoint) & EEG.times <= trialOffsets(thisTimePoint));
                eegData = EEG.data(chanInd, timeInds, thisTrial);
                
                
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
                
                [sustainednessTable] = calcAllFreqPepisodeSustainedness(powerDistribution, frequencies, eegData, EEG.srate, EEG.times, 95, 2);
                sustainednessTable = updateTable(sustainednessTable, subjectID, teleporter, chanNames{thisChan}, thisTrial, timePointNames{thisTimePoint}, spaceType, timeType);
                
                if thisDepth == 1 && thisChan == 1 && thisTrial == 1
                    allSustainednessTable = sustainednessTable;
                else
                    allSustainednessTable = [allSustainednessTable; sustainednessTable];
                end
                
            end % thisTimePoint
            
        end % thisTrial
        
        fprintf('\n\n');
        
        
        
    end % thisChan
    
end % thisDepth

% write the table out to file
fprintf('\n\nWriting data to file...');

writetable(allSustainednessTable, [analysisDir 'csv/' saveFile]);

end

function sustainednessTable = updateTable(sustainednessTable, subjectID, teleporter, electrode, trial, timepoint, spaceType, timeType)
    
    nrow = height(sustainednessTable);
    infoTable = table(repmat({subjectID}, nrow, 1), ...
        repmat({teleporter}, nrow, 1), ...
        repmat({electrode}, nrow, 1), ...
        repmat({trial}, nrow, 1), ...
        repmat({timepoint}, nrow, 1), ...
        repmat({spaceType}, nrow, 1), ...
        repmat({timeType}, nrow, 1));
    infoTable.Properties.VariableNames = {'SubjectID', 'Teleporter', 'Electrode', 'Trial', 'TimePoint', 'SpaceType', 'TimeType'};
    sustainednessTable = horzcat(infoTable, sustainednessTable);

end
