function Teleporter_Epoched_Power(...
    subjectID, ...
    subjectDir, ...
    teleporter, ...
    chanList, ...
    cleanedEpochedPrefix, ...
    cleanedEpochedSuffix, ...
    saveStem, ...
    analysisDir, ...
    saveFile, ...
    calcPower, ...
    timePointNames, ...
    timesNT, ...
    timesFT, ...
    frequencies)
% function Teleporter_Epoched_Power(...
%     subjectID, ...
%     subjectDir, ...
%     teleporter, ...
%     chanList, ...
%     cleanedEpochedPrefix, ...
%     cleanedEpochedSuffix, ...
%     saveStem, ...
%     saveFile, ...
%     calcPower, ...
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
% cleanedEpochedPrefix = {[subjectDir 'Epoched Data/' subjectID '_' teleporter '_epoched_']};
% cleanedEpochedSuffix = {'_noSpikes_noWaves.set'};
% 
% % Specify path to save the power calculations to
% saveStem = [subjectDir 'Mat Files/Power/' subjectID '_' teleporter '_epoched_power_'];
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
% % if you've already calculated power and saved the .mat files, set this to
% % 0, otherwise 1
% calcPower = 0;
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
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/PepisodeCode/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/arne_code/'));

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

% Get depth names from chanList
depthNames = unique(cellfun(@(s) s(1:3), chanList, 'UniformOutput', false));

%% Calculate power for each epoch of each trial
if calcPower
    for thisDepth = 1:length(depthNames)
        
        % Load EEG data
        eegPath = [cleanedEpochedPrefix{1} depthNames{thisDepth} cleanedEpochedSuffix{1}];
        EEG = pop_loadset(eegPath);
        
        %% Create tables of epoch onsets and offsets
        
        % Keep only the 2nd value of each trial type since this indicates the time
        % type (1 = NT, 2 = FT)
        trialTypeList = {EEG.event.type}';
        trialTypeList = cellstr(cellfun(@(s) s(2), trialTypeList));
        
        % Make a table of trial type list
        trialTypeTable = table(trialTypeList);
        
        % Make tables of onsets and offsets
        onsetTimeTable  = table({'1'; '2'}, [timesNT(:, 1)'; timesFT(:, 1)'], 'VariableNames', {'trialTypeList', 'onsetTimes'});
        offsetTimeTable = table({'1'; '2'}, [timesNT(:, 2)'; timesFT(:, 2)'], 'VariableNames', {'trialTypeList', 'offsetTimes'});
        
        % Join the tables together
        trialTimesTable = join(trialTypeTable, onsetTimeTable);
        trialTimesTable = join(trialTimesTable, offsetTimeTable);
        
        
        % Find the channels on this depth electrode
        chanDepth = char(chanList);
        chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
        chanDepth = cellstr(chanDepth);
        chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
        chanNames = chanList(chanInd);
        
        
        %% Calculate power
        calcEpochedPowerLKV(EEG, chanNames, {EEG.event.type}, trialTimesTable{:, 2}, trialTimesTable{:, 3}, saveStem, frequencies, 6, 3);
        
        
    end % thisDepth
end


%% Extract power values for our epochs of interest

% Initialize the cell array to hold all of our pepisode values
powerSummary      = cell(1,10);
powerSummary(1,:) = {'SubjectID','Teleporter','Electrode','TrialNumber','TrialSpaceType','TrialTimeType','TrialType','TimePoint','Frequency','Power'};
thisRow = 2;

for thisDepth = 1:length(depthNames)
    
    fprintf(['\n\nWorking on electrode #' num2str(thisDepth) ' of ' num2str(length(depthNames)) '\n']);
    
    % Find the channels on this depth electrode
    chanDepth = char(chanList);
    chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
    chanDepth = cellstr(chanDepth);
    chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
    chanNames = chanList(chanInd);
    
    % Loop through channels on this depth electrode
    for thisChan = 1:length(chanNames)
        
        % Path to the power vector we calculated in the previous step
        powerVectorFile = [saveStem chanNames{thisChan} '.mat'];
        
        % Load the power vector file
        load(powerVectorFile);
        
        % Loop through trials
        for thisTrial = 1:size(trialTypes, 2)
            
            % Extract the trial type for this trial
            thisLabel = trialTypes{thisTrial};
            
            % Determine whether it's NS or FS
            if strcmpi('1', thisLabel(1)) == 1
                thisSpaceType = 'NS';
            elseif strcmpi('2', thisLabel(1)) == 1
                thisSpaceType = 'FS';
            else
                error('Unknown trial type')
            end
            
            % Determine whether it's NT or FT
            if strcmpi('1', thisLabel(2)) == 1
                thisTimeType = 'NT';
            elseif strcmpi('2', thisLabel(2)) == 1
                thisTimeType = 'FT';
            else
                error('Unknown trial type')
            end
            
            % Combine them together to make the spatiotemporal type
            thisType = [thisSpaceType thisTimeType];
            
            % Loop through time points
            for thisTimePoint = 1:length(timePointNames)
                
                powerVec = allPowerData{thisTrial, thisTimePoint};
                
                % Take the mean across time
                meanPower = nanmean(powerVec, 2);
                
                % Add the values to the summary cell array
                for thisFreq = 1:length(frequencies)
                    
                    powerSummary(thisRow,:) = {subjectID, teleporter, chanNames{thisChan}, thisTrial, thisSpaceType, thisTimeType, thisType, timePointNames{thisTimePoint}, frequencies(thisFreq), meanPower(thisFreq)};
                    
                    thisRow = thisRow + 1;
                    
                end % thisFreq
                
            end % thisTimePoint
            
        end % thisTrial
        
    end % thisChan
    
    
end % thisDepth

% Write out the summary cell array to file
fprintf('\n\nWriting output to file...\n');
dlmcell([analysisDir 'csv/' saveFile '.csv'], powerSummary, 'delimiter', ',');
save([analysisDir 'mat/' saveFile '.mat'], 'powerSummary');
