function [analyses, params] = phaseLockingAnalysis(epochedEEGPath, goodChans, frequencies, waveletCycles, rayleighCycles, rayleighThresh, visualize)
% phaseLockingAnalysis: test whether significant phase locking is observed
% in the 500 ms after teleporter entry or exit
% >> analyses = phaseLockingAnalysis(epochedEEGPath, goodChans, frequencies, waveletCycles, rayleighCycles, rayleighThresh, visualize)
%
% Inputs:
%   epochedEEGPath: path to the epoched EEG dataset (.set file)
%   frequencies: vector of frequencies in Hz
%   goodChans: cell array of channel names for analysis
%   waveletCycles: number of cycles to use for wavelet convolution when
%       estimating phase (recommended = 6)
%   rayleighCycles: number of cycles of data required for phase analysis;
%       for example, 500 ms of data with rayleighCycles = 2 means that
%       frequencies < 4 Hz will be excluded
%   rayleighThresh: pvalue threshold for determining significance
%       (conservative is good since we'll be doing many statistical tests)
%   visualize: either 0 or 1; if 1, this will produce plots for each
%       electrode for NT & FT showing a time-frequency plot of intertrial
%       phase consistency; warning, takes a long time to produce plots
%
% Output:
%   analyses: data structure containing the following fields
%       AnalysisName: e.g., 'NT_Entry'
%       InputData: matrix of input phase data
%       TimeInterval: start and end times that were tested
%       GoodFreq: frequencies that were included in the analysis
%       BadFreq: frequencies that were excluded from analysis
%       AnalyzedPhaseData: subset of InputData, restricted to valid time
%           points and frequencies
%       RayleighP: matrix of P values from the Rayleigh test
%       SignifP: binary vector indicating whether RayleighP exceeded
%           rayleighThresh
%
% Lindsay Vass
% 15 September 2015

% FOR TESTING PURPOSES
% epochedEEGPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Epoched Data/UCDMC15_TeleporterA_epoched_LAD_noSpikes_noWaves.set';
% goodChans      = {'LAD3', 'LAD4', 'LAD5'};
% frequencies    = logspace(log(1)/log(10),log(181)/log(10),31);
% waveletCycles  = 6;
% visualize      = 0; % produce plots of inter-trial phase clustering for each electrode? warning: takes a LONG time
% rayleighCycles = 2; % number of cycles for which you must have significant rayleigh scores
% rayleighThresh = 0.001; % p value must be lower than this to be deemed significant


% calculate instantaneous phase at each time point for each frequency
[ntPhase, ftPhase] = calcInstPhase(epochedEEGPath, frequencies, waveletCycles);

% keep only good channels
EEG = pop_loadset(epochedEEGPath);
chanList = {EEG.chanlocs.labels};
chanInds = nan(size(goodChans));
for i = 1:length(goodChans)
    chanInds(i) = find(strcmpi(goodChans{i}, chanList) == 1);
end

ntPhase = ntPhase(chanInds, :, :, :);
ftPhase = ftPhase(chanInds, :, :, :);
%% make analysis structure
analyses = struct('AnalysisName', NaN, 'Channels', NaN, 'InputData', NaN, 'TimeInterval', NaN, 'GoodFreq', NaN, 'BadFreq', NaN, 'AnalyzedPhaseData', NaN, 'RayleighP', NaN, 'SignifP', NaN);

analyses(1).AnalysisName = 'NT_Entry';
analyses(1).InputData    = ntPhase;
analyses(1).TimeInterval = [0 500];
analyses(1).Channels     = goodChans;

analyses(2).AnalysisName = 'NT_Exit';
analyses(2).InputData    = ntPhase;
analyses(2).TimeInterval = [1830 2331];
analyses(2).Channels     = goodChans;

analyses(3).AnalysisName = 'FT_Entry';
analyses(3).InputData    = ftPhase;
analyses(3).TimeInterval = [0 500];
analyses(3).Channels     = goodChans;

analyses(4).AnalysisName = 'FT_Exit';
analyses(4).InputData    = ftPhase;
analyses(4).TimeInterval = [2830 3331];
analyses(4).Channels     = goodChans;

%% make params structure
params.frequencies    = frequencies;
params.waveletCycles  = waveletCycles;
params.rayleighCycles = rayleighCycles;
params.rayleighThresh = rayleighThresh;
params.samplingRate   = EEG.srate;

%% run Rayleigh tests

for thisAnalysis = 1:length(analyses)
   [phaseDataSample, goodFreq, badFreq, timeMs] = filterPhaseData(analyses(thisAnalysis).InputData, analyses(thisAnalysis).TimeInterval, EEG.times, frequencies, EEG.srate, rayleighCycles);
   analyses(thisAnalysis).GoodFreq = goodFreq;
   analyses(thisAnalysis).BadFreq  = badFreq;
   analyses(thisAnalysis).AnalyzedPhaseData = phaseDataSample;
   analyses(thisAnalysis).TimeMs = timeMs;
   pData = nan(length(goodChans), size(phaseDataSample, 4), size(phaseDataSample, 2));
   for thisElec = 1:size(phaseDataSample, 1)
       for thisFreq = 1:size(phaseDataSample, 4)
           [p, ~] = multRayleighTest(squeeze(phaseDataSample(thisElec, :, :, thisFreq))');
           pData(thisElec, thisFreq, :) = p;
       end
   end
   analyses(thisAnalysis).RayleighP = pData;
   analyses(thisAnalysis).SignifP   = pData < rayleighThresh;
end


%% OPTIONAL: visualize the inter-trial phase locking
if visualize == 1
    ntITPC = nan(size(ntPhase, 1), size(ntPhase, 2), size(ntPhase, 4));
    ftITPC = nan(size(ftPhase, 1), size(ftPhase, 2), size(ftPhase, 4));
    for thisElec = 1:size(ntPhase, 1)
        ntData = squeeze(ntPhase(thisElec, :, :, :));
        ftData = squeeze(ftPhase(thisElec, :, :, :));
        
        ntITPC(thisElec, :, :) = squeeze(abs(mean(exp(1i*ntData), 2)));
        ntITPC = squeeze(ntITPC(thisElec, :, :))';
        ftITPC(thisElec, :, :) = squeeze(abs(mean(exp(1i*ftData), 2)));
        ftITPC = squeeze(ftITPC(thisElec, :, :))';
        
        ntSigP = zeros(size(ntITPC));
        [~, ~, freqInds] = intersect(analyses(1).GoodFreq, frequencies);
        [~, ~, timeInds] = intersect(analyses(1).TimeMs, EEG.times);
        tmpNtSigP = squeeze(analyses(1).SignifP(thisElec, :, :));
        ntSigP(min(freqInds):max(freqInds),min(timeInds):max(timeInds)) = tmpNtSigP;
        ntBounds = bwboundaries(ntSigP);
        
        figure;
        contourf(EEG.times, frequencies, ntITPC, 40, 'linecolor', 'none');
        set(gca, 'ytick', round(logspace(log10(frequencies(1)), log10(frequencies(end)), 10) * 100) / 100, ...
            'yscale', 'log', ...
            'clim', [0 1], ...
            'xlim', [-1000 3000])
        title(['NT Phase Clustering for ' EEG.chanlocs(thisElec).labels])
        ylims = get(gca, 'ylim');
        hold on
        plot([0 0], ylims, 'w--', 'LineWidth', 2)
        plot([1830 1830], ylims, 'w--', 'LineWidth', 2)
        for b = 1:length(ntBounds)
            boundary = ntBounds{b};
            boundary(:,1) = frequencies(boundary(:,1));
            boundary(:,2) = EEG.times(boundary(:,2));
            plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 5);
        end
        hold off
        colorbar
        
        ftSigP = zeros(size(ftITPC));
        [~, ~, freqInds] = intersect(analyses(3).GoodFreq, frequencies);
        [~, ~, timeInds] = intersect(analyses(3).TimeMs, EEG.times);
        tmpFtSigP = squeeze(analyses(3).SignifP(thisElec, :, :));
        ftSigP(min(freqInds):max(freqInds),min(timeInds):max(timeInds)) = tmpFtSigP;
        ftBounds = bwboundaries(ftSigP);
        
        figure;
        contourf(EEG.times, frequencies, ftITPC, 40, 'linecolor', 'none');
        set(gca, 'ytick', round(logspace(log10(frequencies(1)), log10(frequencies(end)), 10) * 100) / 100, ...
            'yscale', 'log', ...
            'clim', [0 1], ...
            'xlim', [-1000 5000])
        title(['FT Phase Clustering for ' EEG.chanlocs(thisElec).labels])
        ylims = get(gca, 'ylim');
        hold on
        plot([0 0], ylims, 'w--', 'LineWidth', 2)
        plot([2830 2830], ylims, 'w--', 'LineWidth', 2)
        for b = 1:length(ftBounds)
            boundary = ftBounds{b};
            boundary(:,1) = frequencies(boundary(:,1));
            boundary(:,2) = EEG.times(boundary(:,2));
            plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 5);
        end
        hold off
        colorbar
    end
end