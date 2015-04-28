% This script will segment the free exploration period into 10 second
% epochs, allow for the removal of bad epochs, and then calculate the mean
% power at each frequency across all epochs, to be used as a baseline in
% future calculations.
%
% Lindsay Vass 27 April 2015

clear all;close all;clc;

%% Load data and set parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
eeglab;

% paths
subject_dir      = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/';
unepochedEEG     = [subject_dir 'Raw Data/UCDMC14_teleporterB_unepoched.set'];
unityFile        = [subject_dir 'Behavioral Data/TeleporterB/s2_patientTeleporterData 2.txt'];
saveFile         = [subject_dir 'Mat Files/UCDMC14_TeleporterB_FreeExplore_baseline.mat']; % mat file with mean power at each electrode/frequency
epochEEGSaveFile = [subject_dir 'Epoched Data/UCDMC14_TeleporterB_FreeExplore_baseline_epoched.mat']; % EEG file of the epoched free exploration data
timeSyncFile     = [subject_dir 'Mat Files/UCDMC14_TeleporterB_time_sync.mat'];

% length of epochs to subdivide data into
epochSec = 10; % in seconds

% load EEG data
EEG = pop_loadset(unepochedEEG);
eeglab redraw;

samplerate  = 512;
widthcycles = 6;
F = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011
wavelet_cycles = 6;

%% Get start times for each mini-epoch

% parse the unity file to get the start and end times of free exploration
% in ticks
fid  = fopen(unityFile);
%  data = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n'); % use this version for unity output that we fixed in Matlab
data =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1); % use this version for raw unity output
fclose(fid);
systemTime = data{3};

% system time in ticks
startTicks = systemTime(1);
endTicks   = systemTime(end);

% convert from ticks to EEG bins
load(timeSyncFile);
startBin = round(startTicks * time_sync_regression(1) + time_sync_regression(2));
endBin   = round(endTicks * time_sync_regression(1) + time_sync_regression(2));

% get epoch length in bins
epochBins = epochSec * EEG.srate;

% how many epochs we'll have
totalEpochs = floor((endBin - startBin) / epochBins);

% create a list of start times for each epoch
epochInds    = nan(totalEpochs,1);
epochInds(1) = startBin;

for thisEpoch = 2:totalEpochs
    
    epochInds(thisEpoch) = epochInds(thisEpoch - 1) + epochBins;
    
end

%% epoch the EEG data

% Insert events into EEG
for n = 1:totalEpochs
    
    EEG.event(n).latency = epochInds(n);
    EEG.event(n).type = '0';
    
end

% Create epochs
[EEG] = pop_epoch(EEG,{},[0 epochSec]);
eeglab redraw;
pop_eegplot(EEG,1,1,1);

% view the data and reject bad epochs
fprintf(' \n\n Remove bad epochs in EEGlab.\n\n Click on bad epochs to highlight them for rejection.\n To unselect an epoch, click it again.\n When done, press REJECT in bottom right.\n')

keyboard;

% Save the epoched data
EEG = pop_saveset(EEG, 'filename', epochEEGSaveFile);

% make sure that the saved data set has fewer epochs than we started with
% and warn us if not
if size(EEG.data, 3) == totalEpochs
    
    answer = questdlg('WARNING! Data size has the same number of epochs we started with. If you rejected any epochs, this means the data did not save correctly. Stop now?','WARNING');
    if strcmpi(answer, 'Yes') == 1
        return
    end
end


%% Calculate mean power in each frequency for each electrode

% other wavelet parameters
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet      = length(time);
n_data         = size(EEG.data,2);
n_convolution  = n_wavelet+n_data-1;
n_conv_pow2    = pow2(nextpow2(n_convolution));

% initalize storage array
meanPower = nan(size(EEG.data,1), size(EEG.data,3), length(F)); % electrodes x epochs x frequencies

for thisChan = 1:size(EEG.data,1) % loop through electrodes
    
    for thisEpoch = 1:size(EEG.data,3) % loop through epochs
    
        % tell us which epoch we're on
        fprintf(['\nComputing mean power for channel ' num2str(thisChan) ' of ' num2str(size(EEG.data,1)) ', Epoch ' num2str(thisEpoch) ' of ' num2str(size(EEG.data,3)) '\n']);
        
        eegdata = squeeze(EEG.data(thisChan,:,thisEpoch));
        
        % FFT of data (note: this doesn't change on frequency iteration)
        fft_data = fft(eegdata,n_conv_pow2);
        
        % initialize output time-frequency data
        tf_data = zeros(length(F),size(eegdata,2));
        
        for fi=1:length(F)
            
            % create wavelet and get its FFT
            wavelet = (pi*F(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*F(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*F(fi)))^2))/F(fi);
            fft_wavelet = fft(wavelet,n_conv_pow2);
            
            % run convolution
            convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
            convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
            
            % identify and exclude outliers, defined as power greater than mean + 3SD
%             temp_data = abs(convolution_result_fft).^2;
            %         orig_data(fi,:) = temp_data;
            %         SD3 = mean(temp_data) + 3*std(temp_data);
            %         ind = find(temp_data >= SD3);
            %         temp_data(ind) = nan;
            
            tf_data(fi,:) = abs(convolution_result_fft).^2;
            
        end
        
        % take mean over time
        mean_tf = nanmean(tf_data,2);
        
        % get median value
        %     med_tf = nanmedian(tf_data,2);
        
        % put data in summary array
        %     medianPower(thisChan,:) = med_tf;
        meanPower(thisChan,thisEpoch,:) = mean_tf;
    end
end

% take the mean across epochs
meanPower = squeeze(mean(meanPower, 2));

% save the results
chanNames = {EEG.chanlocs.labels};
save(saveFile,'meanPower', 'chanNames');