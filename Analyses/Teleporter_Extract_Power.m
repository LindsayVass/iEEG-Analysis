% This script will calculate the power at each frequency for each time
% point of each epoch. It will save this output in a mat file with a
% structure: Timepoint (Pre/Tele/Post) x Condition (NSNT/NSFT/FSNT/FSFT) x
% Electrode x Frequency.
%
% Lindsay Vass 6 May 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% set parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
eeglab;

subjectDir     = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/';
epochedEEGFile = [subjectDir 'Epoched Data/UCDMC14_TeleporterA_epoched.set'];
saveFile       = [subjectDir 'Mat Files/UCDMC14_TeleporterA_power.mat'];

% time periods of interest
preTele = [-1000 0];
teleNT  = [1 1830];
teleFT  = [1 2830];
postNT  = [1831 2830];
postFT  = [2831 3830];

% wavelet parameters (31 log-spaced frequencies, as in Watrous 2011)
waveletcycles = 6;
min_freq = 1;
max_freq = 181;
num_frex = 31;

% channels to analyze
chanList = {'LAD1' 'LHD1' 'RAD1' 'RHD1' 'RHD2'};

% trial types
trialTypeList = {'NSNT' 'NSFT' 'FSNT' 'FSFT'};
trialCodeList = {'11' '12' '21' '22'}; % labels in the EEG events

%% Load data
EEG = pop_loadset(epochedEEGFile);


%% Convert time periods of interest to indices
preTeleInd = find(EEG.times >= preTele(1) & EEG.times <= preTele(2));
teleNTInd  = find(EEG.times >= teleNT(1) & EEG.times <= teleNT(2));
teleFTInd  = find(EEG.times >= teleFT(1) & EEG.times <= teleFT(2));
postNTInd  = find(EEG.times >= postNT(1) & EEG.times <= postNT(2));
postFTInd  = find(EEG.times >= postFT(1) & EEG.times <= postFT(2));

%% Get rid of channels we won't analyze
chans = [];
for thisChan = 1:size(EEG.data,1)
    thisChanName = {EEG.chanlocs(thisChan).labels};
    goodChanInd = strcmpi(thisChanName,chanList);
    if sum(goodChanInd(:)) > 0
        chans(end+1) = thisChan;
    end
end


data = EEG.data(chans,:,:); % keep only channels of interest

% Index epochs by trial type
trialTypes = {EEG.event.type};
trialInds  = cell(length(trialTypeList),1);

for thisType = 1:length(trialTypeList)
    indHolder = [];
    indHolder = find(strcmpi(trialTypes, trialCodeList{thisType}) == 1);
    trialInds{thisType} = indHolder;
end

%% Cohen code

% other wavelet parameters
frequencies = logspace(log(min_freq)/log(10),log(max_freq)/log(10),num_frex);
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet      = length(time);
n_data         = EEG.pnts;
n_convolution  = n_wavelet+n_data-1;
n_conv_pow2    = pow2(nextpow2(n_convolution));

%% Calculate power at each time point for each electrode
allEpochDataDimNames = {'Electrodes x Epochs x Frequencies'};
for thisType = 1:length(trialTypeList)
    
    % get epoch inds for this trial type
    epochInds = trialInds{thisType};
    
    % initialize storage arrays
    allPreData  = nan(length(chans), length(epochInds), length(frequencies));
    allTeleData = allPreData;
    allPostData = allPreData;
    
    for thisChan = 1:length(chans)
        
        % initialize array to hold time-frequency data for all epochs of
        % this condition. epochs x frequencies x time
        allEpochsData = nan(size(trialInds{thisType,1},2), length(frequencies), n_data);
        
        for thisEpoch = 1:size(trialInds{thisType,1},2)
            
            % FFT of data (note: this doesn't change on frequency iteration)
            fft_data = fft(squeeze(data(thisChan,:,epochInds(thisEpoch))),n_conv_pow2);
            
            % initialize output time-frequency data
            tf_data = zeros(length(frequencies),EEG.pnts);
            
            for thisFreq=1:length(frequencies)
                
                % create wavelet and get its FFT
                wavelet = (pi*frequencies(thisFreq)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(thisFreq).*time) .* exp(-time.^2./(2*( waveletcycles /(2*pi*frequencies(thisFreq)))^2))/frequencies(thisFreq);
                fft_wavelet = fft(wavelet,n_conv_pow2);
                
                % run convolution
                convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
                convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
                convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
                
                % put power data into time-frequency matrix
                tf_data(thisFreq,:) = abs(convolution_result_fft).^2;
                
            end % thisFreq
            
            % put data in summary array
            allEpochsData(thisEpoch,:,:) = tf_data(:,:);
            
        end % thisEpoch
        
        % Separate the data into the three time points (Pre / Tele / Post)
        preData = allEpochsData(:, :, preTeleInd);
        
        if ~isempty(strfind(trialTypeList{thisType}, 'NT'))
            teleData = allEpochsData(:, :, teleNTInd);
            postData = allEpochsData(:, :, postNTInd);
        else
            teleData = allEpochsData(:, :, teleFTInd);
            postData = allEpochsData(:, :, postFTInd);
        end
        
        % Take the mean across timepoints
        meanPreData  = squeeze(mean(preData, 3));
        meanTeleData = squeeze(mean(teleData, 3));
        meanPostData = squeeze(mean(postData, 3));
        
        % add to summary array
        allPreData(thisChan, :, :)  = meanPreData;
        allTeleData(thisChan, :, :) = meanTeleData;
        allPostData(thisChan, :, :) = meanPostData;
        
    end % thisChan
    
    tempdata = {allPreData; allTeleData; allPostData};
    allEpochData(thisType) = {tempdata};
    
end % thisType

% save data for all epochs
save(saveFile, 'allEpochData','allEpochDataDimNames','frequencies','trialTypeList','chanList');

