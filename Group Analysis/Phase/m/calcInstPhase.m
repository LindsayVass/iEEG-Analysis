function [ntPhase, ftPhase] = calcInstPhase(epochedEEGPath, frequencies, waveletCycles)
% calcInstPhase: return the phase at each frequency for each timepoint of EEG
% data after convolution with Morlet wavelets; NT and FT trial output is
% stored separately
% >> [ntPhase, ftPhase] = calcInstPhase(epochedEEGPath, frequencies, waveletCycles)
%
% Inputs:
%   epochedEEGPath: path to the epoched EEG dataset (.set file)
%   frequencies: vector of frequencies in Hz at which to extract phase
%   waveletCycles: number of cycles of Morlet wavelets to use for
%       convolution (6 is recommended)
%
% Outputs:
%   ntPhase: phase values for NT trials stored in an electrodes x
%       timepoints x trials x frequencies matrix
%   ftPhase: phase values for FT trials stored in an electrodes x
%       timepoints x trials x frequencies matrix
%
% Lindsay Vass
% 15 September 2015


% FOR TESTING PURPOSES
% epochedEEGPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Epoched Data/UCDMC15_TeleporterA_epoched_LAD_noSpikes_noWaves.set';
% frequencies    = logspace(log(1)/log(10),log(181)/log(10),31);
% waveletCycles  = 6;

% load EEG data
EEG = pop_loadset(epochedEEGPath);

% organize trials by time type
trialList = {EEG.event.type};
trialList = cell2mat(trialList');
trialList = trialList(:, 2); % keep only time label

ntTrials = find(trialList == '1');
ftTrials = find(trialList == '2');

ntPhase = nan(size(EEG.data, 1), size(EEG.data, 2), length(ntTrials), length(frequencies));
ftPhase = nan(size(EEG.data, 1), size(EEG.data, 2), length(ftTrials), length(frequencies));

%% NT trials
% compute instantaneous phase at each time point / frequency
for thisElec = 1:size(EEG.data, 1)
    for thisNT = 1:length(ntTrials)
        eegData  = squeeze(EEG.data(thisElec, :, ntTrials(thisNT)));
        for thisFreq = 1:length(frequencies)
            [~, pha] = morletPowerPhase(eegData, frequencies(thisFreq), EEG.srate, waveletCycles);
            ntPhase(thisElec, :, thisNT, thisFreq) = pha;
        end
    end
end

%% FT trials
% compute instantaneous phase at each time point / frequency
for thisElec = 1:size(EEG.data, 1)
    for thisFT = 1:length(ftTrials)
        eegData  = squeeze(EEG.data(thisElec, :, ftTrials(thisFT)));
        for thisFreq = 1:length(frequencies)
            [~, pha] = morletPowerPhase(eegData, frequencies(thisFreq), EEG.srate, waveletCycles);
            ftPhase(thisElec, :, thisFT, thisFreq) = pha;
        end
    end
end

end

