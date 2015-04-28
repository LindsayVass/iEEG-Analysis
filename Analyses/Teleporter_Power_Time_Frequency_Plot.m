% This script will produce a mean time-frequency plot for each condition for
% each electrode, using the FreeExplore baseline.
%
% Lindsay Vass 28 April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% set parameters
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b'));
eeglab;

subjectDir     = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/';
epochedEEGFile = [subjectDir 'Epoched Data/UCDMC15_TeleporterA_epoched.set'];
baselineFile   = [subjectDir 'Mat Files/UCDMC15_TeleporterA_FreeExplore_baseline.mat'];
saveStem       = [subjectDir 'Figures/UCDMC15_TeleporterA_'];

% time range to plot
plottime = [-1000 4000]; % in ms relative to teleporter entry

% wavelet parameters
samplerate  = 512;
widthcycles = 6;
F = logspace(log(1)/log(10),log(181)/log(10),31); % 31 log-spaced frequencies, as in Watrous 2011

% channels to analyze
chanList = {'RAD4' 'RAD5' 'RAD6' 'RHD2' 'RHD3' 'RHD4' 'LAD3' 'LAD4' 'LHD1' 'LHD2' 'LHD3'};

% trial types
trialTypeList = {'NSNT' 'NSFT' 'FSNT' 'FSFT'};
trialCodeList = {'11' '12' '21' '22'}; % labels in the EEG events


%% Prompt user to select frequencies for analysis
[freqSelect, ~] = listdlg('ListString',cellstr(num2str(F'))','PromptString','Select frequencies to plot');

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

% wavelet parameters
min_freq = 1;
max_freq = 181;
num_frex = 31;

% other wavelet parameters
frequencies = F(freqSelect); 
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


%% Generate time-frequency plots for each electrode
default_clim = [-5 5];

for thisChan = 1:length(chans)
    
    % clear user-defined color limits
    clear user_clim;
    
    % initialize array to hold mean time-frequency for each of the trial
    % types
    allTypesData = nan(length(trialTypeList), length(frequencies), n_data);
    
    for thisType = 1:length(trialTypeList)
        
        % initialize array to hold time-frequency data for all epochs of
        % this condition. epochs x frequencies x time
        allEpochsData = nan(size(trialInds{thisType,1},2), length(frequencies), n_data); 
        
        % get epoch inds for this trial type
        epochInds = trialInds{thisType};
        
        for thisEpoch = 1:size(trialInds{thisType,1},2)
            
            % FFT of data (note: this doesn't change on frequency iteration)
            fft_data = fft(squeeze(data(thisChan,:,epochInds(thisEpoch))),n_conv_pow2);
            
            % initialize output time-frequency data
            tf_data = zeros(length(frequencies),EEG.pnts);
            
            for thisFreq=1:length(frequencies)
                
                % create wavelet and get its FFT
                wavelet = (pi*frequencies(thisFreq)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(thisFreq).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(thisFreq)))^2))/frequencies(thisFreq);
                fft_wavelet = fft(wavelet,n_conv_pow2);
                
                % run convolution
                convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
                convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
                convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
                
                % put power data into time-frequency matrix
                tf_data(thisFreq,:) = abs(convolution_result_fft).^2;
                
            end % thisFreq
            
            % baseline correct
            thisChanName = EEG.chanlocs(chans(thisChan)).labels;
            baselineInd = find(strcmpi(thisChanName,chanNames) == 1);
            baseline = squeeze(meanPower(baselineInd,freqSelect))';
            dbconverted = 10*log10(bsxfun(@rdivide,tf_data,baseline));
            
            % put data in summary array
            allEpochsData(thisEpoch,:,:) = dbconverted(:,:);
            
        end % thisEpoch
        
        % take the average across epochs
        meanEpochsData = squeeze(mean(allEpochsData,1));
        
        % put data in summary array
        allTypesData(thisType,:,:) = meanEpochsData;
        
        if thisType == 1 % allow user to interactively change the color limits
            
            % plot the data
            h = figure;
            contourf(EEG.times(plotind),frequencies,meanEpochsData(:,plotind),40,'linecolor','none')
            set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',default_clim)
            title([trialTypeList{thisType} ' ' chanList{thisChan}])
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
            numTries = 1;
            while 1
                c = get(gca,'clim');
                answer = questdlg('Do you want to change the color scale?');
                if strcmpi(answer,'Yes') == 1
                    newC = inputdlg('Specify new color scale: ','Color scale',1,{num2str(c)});
                    set(gca,'clim',str2num(newC{1}));
                    numTries = numTries + 1;
                else
                    if numTries == 1 % if the default clim is good, then use that
                        user_clim = default_clim;
                    else
                        user_clim = str2num(newC{1}); % otherwise save the new user-defined color limit
                    end % if numTries == 1
                    
                    break
                end % adjust clim
                
            end % clim adjustment while loop
            
        else
            % plot the data
            h = figure;
            contourf(EEG.times(plotind),frequencies,meanEpochsData(:,plotind),40,'linecolor','none')
            set(gca,'ytick',round(logspace(log10(frequencies(1)),log10(frequencies(end)),10)*100)/100,'yscale','log','xlim',plottime,'clim',user_clim)
            title([trialTypeList{thisType} ' ' chanList{thisChan}])
            hcb = colorbar;
            colorHandle = get(hcb,'Title');
            colorLabel = 'dB';
            set(colorHandle,'String',colorLabel);
            hold on
            ylims = get(gca,'YLim');
            plot([0 0],ylims,'--w','LineWidth',2);
            plot([1830 1830],ylims,'--w','LineWidth',2);
            hold off
        end % if thisType == 1
        
        % save the plot
        saveas(2,[saveStem trialTypeList{thisType} '_' chanList{thisChan}],'png');
        
    end % thisType
    
    % Plot all conditions for this channel
    h = figure;
    set(h, 'Position', [0 0 2000 500]);
    
    for thisSubplot = 1:length(trialTypeList)
        
        subaxis(1, length(trialTypeList), thisSubplot, 'Spacing', 0.04);
        contourf(EEG.times, frequencies, squeeze(allTypesData(thisSubplot,:,:)), 40, 'linecolor', 'none');
        set(gca, 'ytick', round(logspace(log10(frequencies(1)), log10(frequencies(end)), 10) * 100) / 100, 'yscale', 'log', 'xlim', plottime, 'clim', user_clim);
        title(trialTypeList{thisSubplot});
        
        hold on;
        ylims = get(gca, 'YLim');
        plot([0 0],ylims,'--w','linewidth',2);
        plot([1830 1830],ylims,'--w','linewidth',2);
        hold off;
        
        hcb = colorbar;
        colorHandle = get(hcb, 'Title');
        colorLabel = 'dB';
        set(colorHandle, 'String', colorLabel);
        
    end % thisSubplot
    
    pause(1);
    
    tightfig;
    
    saveas(h,[saveStem chanList{thisChan}],'png');
    
end % thisChan


