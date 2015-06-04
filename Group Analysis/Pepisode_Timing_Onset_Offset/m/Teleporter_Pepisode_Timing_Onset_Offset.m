function Teleporter_Pepisode_Timing_Onset_Offset(...
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
% % Specify path to save the cell array of pepisode values to. Main
% directory is specified by analysisDir:
% analysisDir = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode_Epoched_Equal_Time_Bins/';
%
% Within that directory, the script will save files using the saveFile
% provided:
% saveFile = [subjectID '_' teleporter '_epoched_pepisode_equal_time_bins'];
% Don't include the filetype. It will save a .mat to analysisDir/mat/ and a
% .csv to analysisDir/csv/
%
% time periods of interest in ms relative to teleporter entry
% timePointNames = {'WholeTrial'};
% 
% timesNT = [-3000 4830];
% timesFT = [-3000 5830];
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

if ~exist([analysisDir 'mat/'], 'dir')
    system(['mkdir ' analysisDirNoSpace 'mat/']);
end

if ~exist([analysisDir 'csv/'], 'dir')
    system(['mkdir ' analysisDirNoSpace 'csv/']);
end

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
        
        % Initialize output arrays
        if charRow == 2
            
            lengthNT = round((timesNT(2) - timesNT(1)) / 1000 * 512);
            lengthFT = round((timesFT(2) - timesFT(1)) / 1000 * 512);
            
            dataNT = cell(1, lengthNT + 2);
            dataFT = cell(1, lengthFT + 2);
            
            dataNT(1,1:2) = {'ObservationID', 'ObservationType'};
            dataFT(1,1:2) = {'ObservationID', 'ObservationType'};
            
            dataNT(1, 3:end) = num2cell(EEG.times(EEG.times >= timesNT(1) & EEG.times <= timesNT(2)));
            dataFT(1, 3:end) = num2cell(EEG.times(EEG.times >= timesFT(1) & EEG.times <= timesFT(2)));
            
        
        end
        
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
                
                [episodeMatrix, onsetMatrix, offsetMatrix] = calcEpochedPepisodeOnsetOffset(powerDistribution, frequencies, eegData, EEG.srate, 95, 3);
              
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
                
                for thisFreq = 1:size(episodeMatrix, 1)
                    observationCharacteristics(charRow, :) = { ...
                        observationID, ...
                        subjectID, ...
                        teleporter, ...
                        chanNames{thisChan}, ...
                        thisTrial, ...
                        spaceType, ...
                        timeType, ...
                        [spaceType timeType], ...
                        timePointNames{thisTimePoint}, ...
                        frequencies(thisFreq)
                        };
                    charRow = charRow + 1;
                    
                    if strcmpi('NT', timeType) == 1
                        
                        dataNT(ntRow, 1:2)       = {observationID, 'Episode'};
                        dataNT(ntRow, 3:end)     = num2cell(episodeMatrix(thisFreq, :));
                        dataNT(ntRow + 1, 1:2)   = {observationID, 'Onset'};
                        dataNT(ntRow + 1, 3:end) = num2cell(onsetMatrix(thisFreq, :));
                        dataNT(ntRow + 2, 1:2)   = {observationID, 'Offset'};
                        dataNT(ntRow + 2, 3:end) = num2cell(offsetMatrix(thisFreq, :));
                        
                        ntRow = ntRow + 3;
                        
                    else
                        
                        dataFT(ftRow, 1:2)       = {observationID, 'Episode'};
                        dataFT(ftRow, 3:end)     = num2cell(episodeMatrix(thisFreq, :));
                        dataFT(ftRow + 1, 1:2)   = {observationID, 'Onset'};
                        dataFT(ftRow + 1, 3:end) = num2cell(onsetMatrix(thisFreq, :));
                        dataFT(ftRow + 2, 1:2)   = {observationID, 'Offset'};
                        dataFT(ftRow + 2, 3:end) = num2cell(offsetMatrix(thisFreq, :));
                        
                        ftRow = ftRow + 3;
                        
                    end
                    
                    observationID = observationID + 1;
                    
                end % thisFreq
                
            end % thisTimePoint
            
        end % thisTrial
        
        fprintf('\n\n');
        
    end % thisChan
    
end % thisDepth


% Write out the summary cell array to file
fprintf('\n\nWriting observation characteristics to file...');
celltocsv([analysisDir 'csv/' saveFile '_observation_characteristics.csv'], observationCharacteristics, 1);
fprintf('\n\nWriting NT data to file...');
celltocsv([analysisDir 'csv/' saveFile '_NT_data.csv'], dataNT, 1);
fprintf('\n\nWriting FT data to file...\n');
celltocsv([analysisDir 'csv/' saveFile '_FT_data.csv'], dataFT, 1);
% save([analysisDir 'mat/' saveFile '_observation_characteristics.mat'], 'observationCharacteristics');
% save([analysisDir 'mat/' saveFile '_NT_data.mat'], 'dataNT');
% save([analysisDir 'mat/' saveFile '_FT_data.mat'], 'dataFT');
