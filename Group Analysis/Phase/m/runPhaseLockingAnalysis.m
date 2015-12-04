% This script will run "phaseLockingAnalysis.m" for each depth electrode
% for each subject/session. It will save a .mat file containing a struct
% called "phaseResults" described below.
%
% phaseResults.subjectID: e.g., 'UCDMC13'
% phaseResults.teleporter.name: e.g., 'TeleporterA'
% phaseResults.teleporter.depths
%   .name: e.g., 'LAD'
%   .chanList: channels to analyze, e.g., 'LAD1, LAD2'
%   .epochedEEGPath: path to the epoched EEG data
%   .results
%       .AnalysisName: denotes time and teleporter condition, e.g.,
%           'NT_Entry' or 'FT_Exit'
%       .Channels: channels that were analyzed
%       .InputData: electrodes x timepoints x trials x frequencies matrix
%           of phase data
%       .TimeInterval: start and end times of analysis in ms, relative to
%           teleporter entry
%       .GoodFreq: list of frequencies valid for analysis
%       .BadFreq: list of frequencies excluded from analysis because they
%           are too slow to complete x cycles in 500 ms, where x is
%           specified by rayleighCycles below
%       .AnalyzedPhaseData: electrodes x timepoints x trials x frequencies
%           matrix of phase data (this is a subset of InputData, restricted
%           in time and frequency)
%       .RayleighP: electrodes x frequencies x timepoints matrix of p
%           values associated with the Rayleigh statistic at each
%           timepoint; significant P values indicate that the distribution
%           of phases across trials was not uniform
%       .SignifP: electrodes x frequencies x timepoints binary matrix
%           indicating whether each P value in RayleighP exceeded the
%           threshold specified by rayleighThresh below
%
% Lindsay Vass
% 16 September 2015

sessionInfoPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo2.mat';

% vector of frequencies in Hz
frequencies = logspace(log(1)/log(10),log(181)/log(10),31);

% number of cycles to use for wavelet convolution when estimating phase 
% (recommended = 6)
waveletCycles = 6;

% number of cycles of data required for phase analysis
% for example, 500 ms of data with rayleighCycles = 2 means that 
% frequencies < 4 Hz will be excluded
rayleighCycles = 2;

% pvalue threshold for determining significance (conservative is good 
% since we'll be doing many statistical tests)
rayleighThresh = 0.01;

% either 0 or 1; if 1, this will produce plots for each electrode for NT &
% FT showing a time-frequency plot of intertrial phase consistency;
% warning, takes a long time to produce plots
visualize = 1;

% where to save the results
savePath = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Phase/mat/phaseResults_rayleighP01.mat';

%% Run analysis

load(sessionInfoPath);
phaseResults = sessionInfo;

for thisSubject = 1:length(sessionInfo)
    for thisTele = 1:length(sessionInfo(thisSubject).teleporter)
        for thisDepth = 1:length(sessionInfo(thisSubject).teleporter(thisTele).depths)
            epochedEEGPath = sessionInfo(thisSubject).teleporter(thisTele).depths(thisDepth).epochedEEGPath;
            goodChans      = sessionInfo(thisSubject).teleporter(thisTele).depths(thisDepth).chanList;
            [analyses, params] = phaseLockingAnalysis(epochedEEGPath, goodChans, frequencies, waveletCycles, rayleighCycles, rayleighThresh, visualize);
            phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results = analyses;
        end
    end
end

save(savePath, 'phaseResults', 'params', '-v7.3');