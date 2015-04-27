% This script will produce a mean time-frequency plot for each condition for
% each electrode, using the FreeExplore stationary baseline.
%
% Lindsay Vass 24 April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% set parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
addpath(genpath('/Users/mcopara/Code/iEEG_code/code/downloaded_code/eeg_toolbox/eeg_toolbox_v1_3_2'));
addpath(genpath('/Users/mcopara/Code/iEEG_code/code/arne_code/matlab_scripts_from_ucla/'));
eeglab;

epochedEEGFile = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Raw Data/UCDMC15_teleporterA_epoched.set';
baselineFile = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Mat Files/UCDMC15_TeleporterA_FreeExplore_stationary_baseline.mat';

savePath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Figures/TeleporterA_';

% time range to plot
plottime = [-1000 4000]; % in ms relative to teleporter entry

% wavelet parameters
samplerate  = 512;
widthcycles = 6;
F = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

% channels to analyze
chanList = {'RAD4' 'RAD5' 'RAD6' 'RHD2' 'RHD3' 'RHD4' 'LAD3' 'LAD4' 'LHD1' 'LHD2' 'LHD3'};


%% Load data

EEG = pop_loadset(epochedEEGFile);
load(baselineFile);

%% Get rid of channels we won't analyze
chans = [];
for thisChan = 1:size(EEG.data,1)
    thisChanName = {EEG.chanlocs(thisChan).labels};
    goodChanInd = strcmpi(thisChanName,chanList);
    if sum(goodChanInd(:)) > 0
        chans(end+1) = thisChan;
    end
end


data = EEG.data;
data = data(chans,:,:); % keep only channels of interest

% Separate out good epochs by trial type
trialTypes = {EEG.event.type};
NSNT = find(strcmpi('11',trialTypes) == 1);
NSFT = find(strcmpi('12',trialTypes) == 1);
FSNT = find(strcmpi('21',trialTypes) == 1);
FSFT = find(strcmpi('22',trialTypes) == 1);


%% Cohen code

% wavelet parameters
min_freq = 1;
max_freq = 181;
num_frex = 31;

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
frequencies = frequencies(1:20); % exclude gamma band (30+ Hz)
time = -1:1/EEG.srate:1;
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet      = length(time);
n_data         = EEG.pnts;
n_convolution  = n_wavelet+n_data-1;
n_conv_pow2    = pow2(nextpow2(n_convolution));
wavelet_cycles = 6;


% define plot period
plotind  = find(EEG.times >= plottime(1) & EEG.times <= plottime(2));



%% Compute power for NSNT
ALL_NSNT = nan(length(chans),20,n_data);
for thisChan = 1:length(chans)
    
    all_NSNT_data = nan(length(NSNT),20,n_data);
    
    for thisEpoch = 1:length(NSNT)
        
        % FFT of data (note: this doesn't change on frequency iteration)
        fft_data = fft(squeeze(EEG.data(thisChan,:,NSNT(thisEpoch))),n_conv_pow2);
        
        % initialize output time-frequency data
        tf_data = zeros(length(frequencies),EEG.pnts);
        
        for fi=1:length(frequencies)
            
            % create wavelet and get its FFT
            wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
            fft_wavelet = fft(wavelet,n_conv_pow2);
            
            % run convolution
            convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
            convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
            
            % put power data into time-frequency matrix
            tf_data(fi,:) = abs(convolution_result_fft).^2;
            
        end
        
        % baseline correct
        thisChanName = EEG.chanlocs(chans(thisChan)).labels;
        baselineInd = find(strcmpi(thisChanName,chanNames) == 1);
        baseline = squeeze(meanPower(baselineInd,1:20))';
        %         baseline_rm = repmat(baseline,[1 size(tf_data,2)]);
        %         prctChange = 100 * (tf_data - baseline_rm)./baseline_rm;
        dbconverted = 10*log10(bsxfun(@rdivide,tf_data,baseline));
        
        
        % put data in summary array
        all_NSNT_data(thisEpoch,:,:) = dbconverted(:,:);
        
    end
    
    
    grand_mean_NSNT_data = squeeze(mean(all_NSNT_data,1));
    
    figure(2);
    contourf(EEG.times(plotind),frequencies,grand_mean_NSNT_data(:,plotind),40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title(['NSNT ' chanList{thisChan}])
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w','LineWidth',2);
    plot([1830 1830],ylims,'--w','LineWidth',2);
    hold off
    
    % allow user to change colorbar settings until it looks good
    while 1
        c = get(gca,'clim');
        answer = questdlg('Do you want to change the color scale?');
        if strcmpi(answer,'Yes') == 1
            newC = inputdlg('Specify new color scale: ','Color scale',1,{num2str(c)});
            set(gca,'clim',str2num(newC{1}));
        else
            break
        end
        
    end
    
    saveas(2,[savePath 'NSNT_Hipp_Power_' chanList{thisChan}],'png');
    
    ALL_NSNT(thisChan,:,:) = grand_mean_NSNT_data(:,:);
end



%% Compute power for NSFT
ALL_NSFT = nan(length(chans),20,n_data);
for thisChan = 1:length(chans)
    all_NSFT_data = nan(length(NSFT),20,n_data);
    
    for thisEpoch = 1:length(NSFT)
        
        % FFT of data (note: this doesn't change on frequency iteration)
        fft_data = fft(squeeze(EEG.data(thisChan,:,NSFT(thisEpoch))),n_conv_pow2);
        
        % initialize output time-frequency data
        tf_data = zeros(length(frequencies),EEG.pnts);
        
        for fi=1:length(frequencies)
            
            % create wavelet and get its FFT
            wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
            fft_wavelet = fft(wavelet,n_conv_pow2);
            
            % run convolution
            convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
            convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
            
            % put power data into time-frequency matrix
            tf_data(fi,:) = abs(convolution_result_fft).^2;
            
        end
        
        % baseline correct
       thisChanName = EEG.chanlocs(chans(thisChan)).labels;
        baselineInd = find(strcmpi(thisChanName,chanNames) == 1);
        baseline = squeeze(meanPower(baselineInd,1:20))';
        %         baseline_rm = repmat(baseline,[1 size(tf_data,2)]);
        %         prctChange = 100 * (tf_data - baseline_rm)./baseline_rm;
        dbconverted = 10*log10(bsxfun(@rdivide,tf_data,baseline));
        
        % put data in summary array
        all_NSFT_data(thisEpoch,:,:) = dbconverted(:,:);
        
        
    end
    
    grand_mean_NSFT_data = squeeze(mean(all_NSFT_data,1));
    
    figure(3)
    contourf(EEG.times(plotind),frequencies,grand_mean_NSFT_data(:,plotind),40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title(['NSFT ' chanList{thisChan}])
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w');
    plot([2830 2830],ylims,'--w');
    hold off
    
    % allow user to change colorbar settings until it looks good
    while 1
        c = get(gca,'clim');
        answer = questdlg('Do you want to change the color scale?');
        if strcmpi(answer,'Yes') == 1
            newC = inputdlg('Specify new color scale: ','Color scale',1,{num2str(c)});
            set(gca,'clim',str2num(newC{1}));
        else
            break
        end
        
    end
    
    saveas(3,[savePath 'NSFT_Hipp_Power_' chanList{thisChan}],'png');
    
    ALL_NSFT(thisChan,:,:) = grand_mean_NSFT_data(:,:);
end




%% Compute power for FSNT
ALL_FSNT = nan(length(chans),20,n_data);

for thisChan = 1:length(chans)
    all_FSNT_data = nan(length(FSNT),20,n_data);
    
    for thisEpoch = 1:length(FSNT)
        
        % FFT of data (note: this doesn't change on frequency iteration)
        fft_data = fft(squeeze(EEG.data(thisChan,:,FSNT(thisEpoch))),n_conv_pow2);
        
        % initialize output time-frequency data
        tf_data = zeros(length(frequencies),EEG.pnts);
        
        for fi=1:length(frequencies)
            
            % create wavelet and get its FFT
            wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
            fft_wavelet = fft(wavelet,n_conv_pow2);
            
            % run convolution
            convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
            convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
            
            % put power data into time-frequency matrix
            tf_data(fi,:) = abs(convolution_result_fft).^2;
            
        end
        
        % baseline correct
        thisChanName = EEG.chanlocs(chans(thisChan)).labels;
        baselineInd = find(strcmpi(thisChanName,chanNames) == 1);
        baseline = squeeze(meanPower(baselineInd,1:20))';
        %         baseline_rm = repmat(baseline,[1 size(tf_data,2)]);
        %         prctChange = 100 * (tf_data - baseline_rm)./baseline_rm;
        dbconverted = 10*log10(bsxfun(@rdivide,tf_data,baseline));
        
        % put data in summary array
        all_FSNT_data(thisEpoch,:,:) = dbconverted(:,:);
        
    end
    
    grand_mean_FSNT_data = squeeze(mean(all_FSNT_data,1));
    
    
    figure(4)
    contourf(EEG.times(plotind),frequencies,grand_mean_FSNT_data(:,plotind),40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title(['FSNT ' chanList{thisChan}])
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w');
    plot([1830 1830],ylims,'--w');
    hold off
    
    % allow user to change colorbar settings until it looks good
    while 1
        c = get(gca,'clim');
        answer = questdlg('Do you want to change the color scale?');
        if strcmpi(answer,'Yes') == 1
            newC = inputdlg('Specify new color scale: ','Color scale',1,{num2str(c)});
            set(gca,'clim',str2num(newC{1}));
        else
            break
        end
        
    end
    saveas(4,[savePath 'FSNT_Hipp_Power_' chanList{thisChan}],'png');
    
    ALL_FSNT(thisChan,:,:) = grand_mean_FSNT_data(:,:);
end





%% Compute power for FSFT
ALL_FSFT = nan(length(chans),20,n_data);
for thisChan = 1:length(chans)
    all_FSFT_data = nan(length(FSFT),20,n_data);
    
    for thisEpoch = 1:length(FSFT)
        
        % FFT of data (note: this doesn't change on frequency iteration)
        fft_data = fft(squeeze(EEG.data(thisChan,:,FSFT(thisEpoch))),n_conv_pow2);
        
        % initialize output time-frequency data
        tf_data = zeros(length(frequencies),EEG.pnts);
        
        for fi=1:length(frequencies)
            
            % create wavelet and get its FFT
            wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
            fft_wavelet = fft(wavelet,n_conv_pow2);
            
            % run convolution
            convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
            convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
            
            % put power data into time-frequency matrix
            tf_data(fi,:) = abs(convolution_result_fft).^2;
            
        end
        
        % baseline correct
       thisChanName = EEG.chanlocs(chans(thisChan)).labels;
        baselineInd = find(strcmpi(thisChanName,chanNames) == 1);
        baseline = squeeze(meanPower(baselineInd,1:20))';
        %         baseline_rm = repmat(baseline,[1 size(tf_data,2)]);
        %         prctChange = 100 * (tf_data - baseline_rm)./baseline_rm;
        dbconverted = 10*log10(bsxfun(@rdivide,tf_data,baseline));
        
        % put data in summary array
        all_FSFT_data(thisEpoch,:,:) = dbconverted(:,:);
        
        
    end
    
    grand_mean_FSFT_data = squeeze(mean(all_FSFT_data,1));
    
    figure(5)
    contourf(EEG.times(plotind),frequencies,grand_mean_FSFT_data(:,plotind),40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title(['FSFT ' chanList{thisChan}])
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w');
    plot([2830 2830],ylims,'--w');
    hold off
    
    % allow user to change colorbar settings until it looks good
    while 1
        c = get(gca,'clim');
        answer = questdlg('Do you want to change the color scale?');
        if strcmpi(answer,'Yes') == 1
            newC = inputdlg('Specify new color scale: ','Color scale',1,{num2str(c)});
            set(gca,'clim',str2num(newC{1}));
        else
            break
        end
        
    end
    saveas(5,[savePath 'FSFT_Hipp_Power_' chanList{thisChan}],'png');
    
    ALL_FSFT(thisChan,:,:) = grand_mean_FSFT_data(:,:);
end



%% Plot all 4 conditions for each channel
for thisChan = 1:length(chans)
    
    thisNSNT = squeeze(ALL_NSNT(thisChan,:,:));
    thisNSFT = squeeze(ALL_NSFT(thisChan,:,:));
    thisFSNT = squeeze(ALL_FSNT(thisChan,:,:));
    thisFSFT = squeeze(ALL_FSFT(thisChan,:,:));
    
    h = figure(6);
    set(h,'Position',[0 0 2000 500]);
    
    subaxis(1,4,1,'Spacing',0.04);
    contourf(EEG.times,frequencies,thisNSNT,40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title('NSNT')
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w','linewidth',2);
    plot([1830 1830],ylims,'--w','linewidth',2);
    hold off
    
    subaxis(1,4,2,'Spacing',0.04);
    contourf(EEG.times,frequencies,thisNSFT,40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title('NSFT')
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w','linewidth',2);
    plot([2830 2830],ylims,'--w','linewidth',2);
    hold off
    
    subaxis(1,4,3,'Spacing',0.04);
    contourf(EEG.times,frequencies,thisFSNT,40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title('FSNT')
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w','linewidth',2);
    plot([1830 1830],ylims,'--w','linewidth',2);
    hold off
    
    subaxis(1,4,4,'Spacing',0.04);
    contourf(EEG.times,frequencies,thisFSFT,40,'linecolor','none')
    set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',[-5 5])
    title('FSFT')
    hcb = colorbar;
    colorHandle = get(hcb,'Title');
    colorLabel = 'dB';
    set(colorHandle,'String',colorLabel);
    hold on
    ylims = get(gca,'YLim');
    plot([0 0],ylims,'--w','linewidth',2);
    plot([2830 2830],ylims,'--w','linewidth',2);
    hold off
    
    pause(1);
    
    tightfig;
    
    saveas(h,[savePath 'Hipp_Power_' chanList{thisChan}],'png');
    
end


%% Save data
% save('/Users/lindsay/Documents/MATLAB/iEEG/Group Analysis/TeleporterPower/UCDMC14_Hipp_TF.mat','grand_mean_all_data','grand_mean_FSFT_data','grand_mean_FSNT_data','grand_mean_NSFT_data','grand_mean_NSNT_data');