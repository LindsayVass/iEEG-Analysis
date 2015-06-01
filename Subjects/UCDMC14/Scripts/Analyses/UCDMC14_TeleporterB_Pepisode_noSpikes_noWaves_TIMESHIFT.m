% This script will perform a pepisode analysis, which consists of two
% parts. First, the script will calculate the power in each frequency
% across the entire experiment and use this to build a distribution of
% power values for each frequency. The value at 95% of the distribution is
% used as a power threshold in the second part of the analysis. For each
% epoch, the script will calculate the percent of time spent oscillating at
% a given frequency, where the oscillations must exceed both a power
% threshold (95% of the power distribution) and a duration threshold (# of
% cycles).
%
% In this script, we will evaluate three epochs:
% 1. Pre-teleportation (-1000 : 0 ms)
% 2. Teleportation (0 : 1830 ms for NT or 0 : 2830 ms for FT)
% 3. Post-teleportation (1830 : 1930 ms for NT or 2830 : 2930 ms for FT)
%
% Lindsay Vass 29 April 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% set parameters for analysis

% Subject info
subjectID  = 'UCDMC14';
subjectDir = ['/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/' subjectID '/'];
teleporter = 'TeleporterB';

% Specify file naming conventions for data. There are
% separate EEG files for each depth electrode, so we will specify every
% part of the path except the depth electrode name. Then, we will later
% combine these as [prefix depthName suffix]. We will use a cell array to
% allow for multiple prefixes or suffixes (e.g., different prefixes for
% EDF1 and EDF2)
cleanedUnepochedPrefix = {[subjectDir 'PreProcessing Intermediates/' subjectID '_' teleporter '_unepoched_']};
cleanedUnepochedSuffix = {'_noSpikes_noWaves.set'};

% Do the same for the epoched cleaned data
% cleanedEpochedPrefix   = {[subjectDir 'Epoched Data/' subjectID '_' teleporter '_epoched_']};
% cleanedEpochedSuffix   = {'_noSpikes_noWaves.set'};

% Specify path to save the pepisode calculations to
saveStem = [subjectDir 'Mat Files/Pepisode/' subjectID '_' teleporter '_pepisode_'];

% Specify path to save the cell array of pepisode values to
saveFile = [subjectDir 'Mat Files/Pepisode/Summary/' subjectID '_' teleporter '_pepisode_summary_noSpikes_noWaves_TIMESHIFT'];

% channel names to use
chanList  = {'LAD1' 'LHD1' 'LHD2' 'RAD1' 'RHD1' 'RHD2'};

% thresholds for pepisode 
durationThresh  = 3; % # of cycles
amplitudeThresh = 95; % percent of distribution
    
% if the distribution of power across the recording has already been
% calculated, set this to 1
skipCompute = 1;

% time periods of interest in ms relative to teleporter entry
timePointNames = {'Pre' 'Tele' 'Post'};
preTele = [-1000 0];
teleNT  = [1 1830];
teleFT  = [1 2830];
postNT  = [1831 2830];
postFT  = [2831 3830];

% frequencies to use
frequencies = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

%% set paths and filenames
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/PepisodeCode/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'));
addpath(genpath('/Users/Lindsay/Documents/MATLAB/arne_code/'));

% Make the output directory if it doesn't already exist
if ~exist([subjectDir 'Mat Files/Pepisode/'], 'dir')
    system(['mkdir ' subjectDir 'Mat\ Files/Pepisode/']);
end

if ~exist([subjectDir 'Mat Files/Pepisode/Summary/'], 'dir')
    system(['mkdir ' subjectDir 'Mat\ Files/Pepisode/Summary/']);
end

% Get depth names from chanList
depthNames = unique(cellfun(@(s) s(1:3), chanList, 'UniformOutput', false));

eeglab;


%% Calculate power distributions for pepisode

% In this first step, we will use the cleaned, unepoched EEG files to
% establish a distribution of power values at each frequency of interest.
% This section will return a binary vector  that contains a 0 for each
% point in time without a sustained oscillation at that frequency and a 1
% for each point in time that exceeds both the duration and amplitude
% thresholds (i.e., is a sustained oscillation)

if ~skipCompute
    
    for thisDepth = 1:length(depthNames)
        
        % Find the channels on this depth electrode
        chanDepth = char(chanList);
        chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
        chanDepth = cellstr(chanDepth);
        chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
        chanNames = chanList(chanInd);
        
        % Set up eeg list for this channel
        numEDFs = size(cleanedUnepochedPrefix,2);
        eegList = {};
        for thisEDF = 1:numEDFs
            pathName = [cleanedUnepochedPrefix{thisEDF} depthNames{thisDepth} cleanedUnepochedSuffix{1}];
            eegList{thisEDF} = pathName;
            
        end % thisEDF
        
        % Calculate pepisode
        calcPepisodeLKV(eegList, chanNames, saveStem, frequencies, durationThresh, amplitudeThresh);
        
        
    end % thisDepth
    
end

%% Extract pepisode values for our epochs of interest

% Initialize the cell array to hold all of our pepisode values
pepisodeSummary      = cell(1,11);
pepisodeSummary(1,:) = {'SubjectID','Teleporter','EDF','Electrode','TrialNumber','TrialSpaceType','TrialTimeType','TrialType','TimePoint','Frequency','Pepisode'};
thisRow = 2;

for thisDepth = 1:length(depthNames)
    
    fprintf(['\n\nWorking on electrode #' num2str(thisDepth) ' of ' num2str(length(depthNames)) '\n']);
    
    % Find the channels on this depth electrode
    chanDepth = char(chanList);
    chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
    chanDepth = cellstr(chanDepth);
    chanInd   = strcmpi(depthNames{thisDepth}, chanDepth);
    chanNames = chanList(chanInd);
    
    numEDFs = size(cleanedUnepochedPrefix,2);
    
    % Loop through EEG files
    for thisEDF = 1:numEDFs
        
        % Load cleaned unepoched data
        EEG = pop_loadset([cleanedUnepochedPrefix{thisEDF} depthNames{thisDepth} cleanedUnepochedSuffix{1}]);
        
        % Extract the trial events from the event list (i.e., remove boundary
        % events)
        boundaryInds = strcmpi('boundary', {EEG.event.type});
        trialList    = EEG.event(boundaryInds == 0); % select events that are NOT boundaries
        
        % Loop through channels on this depth electrode
        for thisChan = 1:length(chanNames)
            
            % Path to the binary vector we calculated in the previous step
            binaryVectorFile = [saveStem chanNames{thisChan} '.mat'];
            
            % Loop through trials
            for thisTrial = 1:size(trialList, 2)
                
                % Extract the trial type for this trial
                thisLabel = trialList(thisTrial).type;
                
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
                    
                    % Obtain the binary vector for the timepoint of
                    % interest
                    switch thisTimePoint
                        case 1 % Pre-teleportation
                            
                            binarySelection = getPepisodeLKV(binaryVectorFile, thisEDF, trialList(thisTrial).latency, preTele(2) - preTele(1), preTele(1) + 10000, frequencies);
                        
                        case 2 % Teleportation
                            
                            if strcmpi('NT', thisTimeType) == 1
                                binarySelection = getPepisodeLKV(binaryVectorFile, thisEDF, trialList(thisTrial).latency, teleNT(2) - teleNT(1), teleNT(1) + 10000, frequencies);
                            else
                                binarySelection = getPepisodeLKV(binaryVectorFile, thisEDF, trialList(thisTrial).latency, teleFT(2) - teleFT(1), teleFT(1) + 10000, frequencies);
                            end
                            
                        case 3 % Post-teleportation
                            
                            if strcmpi('NT', thisTimeType) == 1
                                binarySelection = getPepisodeLKV(binaryVectorFile, thisEDF, trialList(thisTrial).latency, postNT(2) - postNT(1), postNT(1) + 10000, frequencies);
                            else
                                binarySelection = getPepisodeLKV(binaryVectorFile, thisEDF, trialList(thisTrial).latency, postFT(2) - postFT(1), postFT(1) + 10000, frequencies);
                            end
                    end % switch thisTimePoint
                    
                    % Take the mean across time
                    meanPepisode = mean(binarySelection,3);
                    
                    % Add the values to the summary cell array
                    for thisFreq = 1:length(frequencies)
                        
                        pepisodeSummary(thisRow,:) = {subjectID, teleporter, thisEDF, chanNames{thisChan}, thisTrial, thisSpaceType, thisTimeType, thisType, timePointNames{thisTimePoint}, frequencies(thisFreq), meanPepisode(thisFreq)};
                        
                        thisRow = thisRow + 1;
                        
                    end % thisFreq
                    
                end % thisTimePoint
                 
            end % thisTrial
            
        end % thisChan
        
    end % thisEDF
    
end % thisDepth

% Write out the summary cell array to file
fprintf('\n\nWriting output to file...\n');
dlmcell([saveFile '.csv'], pepisodeSummary, 'delimiter', ',');
save([saveFile '.mat'], 'pepisodeSummary');










