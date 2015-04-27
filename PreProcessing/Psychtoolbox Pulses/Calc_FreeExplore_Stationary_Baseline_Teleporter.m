% This script will identify the stationary period at the beginning of free
% exploration and calculate the mean power in each frequency in each
% electrode during this period.
%
% Lindsay Vass 27 April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% Load data and set parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
eeglab;

% Paths
subjectDir       = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/';
unityDataPath    = [subjectDir 'Behavioral Data/TeleporterB/s2_freeexplore_patientTeleporterData 2.txt'];
outFile          = [subjectDir 'Mat Files/UCDMC14_TeleporterB_FreeExplore_stationary_baseline.mat'];
timeSyncFile     = [subjectDir 'Mat Files/UCDMC14_TeleporterB_time_sync.mat'];
unepochedEEGFile = [subjectDir 'Raw Data/UCDMC14_teleporterB_unepoched.set'];
EEG = pop_loadset(unepochedEEGFile);
data = EEG.data;

% wavelet parameters
F = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;
wavelet_cycles = 6;


%% identify stationary period

% load data
fid  = fopen(unityDataPath);
txtData = textscan(fid,'%f%f%f%f%f','delimiter',',','Headerlines',1);
fclose(fid);

systemTime = txtData{2};
xPos = txtData{3};
zPos = txtData{4};
yRot = txtData{5};

% start time of free exploration
startTime = systemTime(1);

thisLine = 2;
while 1
    
    % check if x,y,z are same as previous time point
    if xPos(thisLine) == xPos(thisLine - 1) && yRot(thisLine) == yRot(thisLine - 1) && zPos(thisLine) == zPos(thisLine - 1)
        thisLine = thisLine + 1;
    else
        endTime = systemTime(thisLine - 1);
        break
    end
    
end

%% convert startTime and endTime from ticks to EEG bins

load(timeSyncFile);
startBin = round(startTime * time_sync_regression(1) + time_sync_regression(2));
endBin   = round(endTime * time_sync_regression(1) + time_sync_regression(2));

%% Visualize the baseline EEG so we can remove any time periods with noise

% Insert events for baseline start/end
EEG.event(1).latency = startBin;
EEG.event(1).type = '00';
EEG.event(2).latency = endBin;
EEG.event(2).type = '00';
eeglab redraw;
pop_eegplot(EEG,1,1,1);

fprintf('\n\n Check that the baseline period (events 00 mark beginning and end) looks ok.\n\n');

keyboard;

%% Calculate mean power in each frequency across entire recording

% FFT parameters (use next-power-of-2)
n_wavelet      = length(time);
n_data         = endBin - startBin + 1;
n_convolution  = n_wavelet+n_data-1;
n_conv_pow2    = pow2(nextpow2(n_convolution));

chanNames = cell(size(data,1),1);
medianPower = nan(size(data,1),length(F));
for thisChan = 1:size(data,1)
    
    thisChanName = EEG.chanlocs(thisChan).labels;
    chanNames{thisChan} = thisChanName;
    
    fprintf(['computing mean power for channel ' num2str(thisChan) ' of ' num2str(size(data,1)) '\n']);
    eegdata = squeeze(data(thisChan,startBin:endBin));
    
    % FFT of data (note: this doesn't change on frequency iteration)
    fft_data = fft(eegdata,n_conv_pow2);
    
    % initialize output time-frequency data
    tf_data = zeros(length(F),size(eegdata,2));
    orig_data = tf_data;
    
    for fi=1:length(F)
        
        % create wavelet and get its FFT
        wavelet = (pi*F(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*F(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*F(fi)))^2))/F(fi);
        fft_wavelet = fft(wavelet,n_conv_pow2);
        
        % run convolution
        convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
        convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
        convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        % identify and exclude outliers, defined as power greater than mean + 3SD
        temp_data = abs(convolution_result_fft).^2;
        %         orig_data(fi,:) = temp_data;
        %         SD3 = mean(temp_data) + 3*std(temp_data);
        %         ind = find(temp_data >= SD3);
        %         temp_data(ind) = nan;
        
        tf_data(fi,:) = temp_data;
        
    end
    
    % take mean over time
    mean_tf = nanmean(tf_data,2);
    
    meanPower(thisChan,:) = mean_tf;
    
end

% save the results
% save(outFile,'medianPower');
save(outFile,'meanPower','chanNames');