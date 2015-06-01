
subjects = {'UCDMC14' 'UCDMC15'};
versions = {'TeleporterA' 'TeleporterB'};

exp_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/';

outFile = ['/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/Wide_Pepisode_by_condition_' date '.csv'];
rOutFile = ['/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Pepisode/Tall_Pepisode_by_condition_' date '.csv'];

%% load in the data for each subject/version

allData = cell(length(subjects), length(versions));
for thisSubject = 1:length(subjects)
    
    for thisVersion = 1:length(versions)
        
        % navigate to subject's directory
        cd([exp_dir subjects{thisSubject} '/Mat Files/']);
        
        % load the mat file
        load([subjects{thisSubject} '_' versions{thisVersion} '_pepisode_by_condition.mat']);
        
        % put it into the holder array
        allData(thisSubject,thisVersion) = {PepisodeByFreqBand};
        
    end % thisVersion
    
end % thisSubject

%% organize data for R output

rTable = {};
rTable(1,:) = {'Electrode','TrialType','FreqBand','TimePoint','Power'};
thisRow = 1;
freqBandNames = {'Delta' 'Theta' 'Alpha' 'Beta' 'Gamma'};
timePointNames = {'Pre' 'Tele' 'Post'};
varNames = {};
chanNumber = 1;
for thisDataset = 1:size(allData,1) * size(allData,2)
    
    thisData = allData{thisDataset};
    
    for thisChan = 1:size(thisData,1)
        
        for thisType = 1:length(trialTypeList)
            
            for thisFreq = 1:size(PepisodeByFreqBand,3)
                
                for thisTimePoint = 1:size(PepisodeByFreqBand,4)
                    
                    rTable(thisRow, :) = {num2str(chanNumber), trialTypeList{thisType}, freqBandNames{thisFreq}, timePointNames{thisTimePoint}, thisData(thisChan, thisType, thisFreq, thisTimePoint)};
                   
                    thisRow = thisRow + 1;
                    
                end % thisTimePoint
                 
            end % thisFreq
            
        end % thisType
        
        chanNumber = chanNumber + 1;
        
    end % thisChan
    
end % thisDataset

% save the data
dlmcell(rOutFile, rTable, 'delimiter',',');


%% organize the data into a table for output
dataTable = [];
thisRow = 1;
freqBandNames = {'Delta' 'Theta' 'Alpha' 'Beta' 'Gamma'};
timePointNames = {'Pre' 'Tele' 'Post'};
varNames = {};
for thisDataset = 1:size(allData,1) * size(allData,2)
    
    thisCol = 1;
    thisData = allData{thisDataset};
    
    for thisChan = 1:size(thisData,1)
        
        for thisType = 1:length(trialTypeList)
            
            for thisFreq = 1:size(PepisodeByFreqBand,3)
                
                for thisTimePoint = 1:size(PepisodeByFreqBand,4)
                    
                    dataTable(thisRow, thisCol) = thisData(thisChan, thisType, thisFreq, thisTimePoint);
                    
                    
                    if thisDataset == 1
                        varNames(1,thisCol) = {[trialTypeList{thisType} '_' freqBandNames{thisFreq} '_' timePointNames{thisTimePoint}]};
                    end
                    
                    thisCol = thisCol + 1;
                    
                end % thisTimePoint
                 
            end % thisFreq
            
        end % thisType
        
        thisRow = thisRow + 1;
        thisCol = 1;
        
    end % thisChan
    
end % thisDataset

% make cell array
csvOutput = cell(size(dataTable,1) + 1, size(dataTable,2));
csvOutput(1,:) = varNames;
csvOutput(2:end,:) = num2cell(dataTable);

% save output
dlmcell(outFile,csvOutput,'delimiter',',');