function calcEpochedPowerLKV(EEG, chanNames, trialTypes, epochOnsets, epochOffsets, saveStem, frequencies, width, bufferCycles)
% calcPowerLKV: Calculate power on the clean segments of the EEG.
% This function will loop through the channels in "chanNames" and calculate
% the z-score of the power of the EEG data at each frequency.
% The function will first loop through all EEG files and calculate the
% power at each time point. It will then log-transform the power and
% normalize it by z-scoring within each electrode and frequency with
% respect to the mean and standard deviation. It will save these values to
% an output file where they can be accessed later using getPowerLKV.
%
% INPUTS:
%   EEG: EEG structure from EEGLAB
%
%   chanNames: cell array of channel names to analyze for this EEG
%
%   trialTypes: cell array of trial type labels; this isn't actually used,
%       just saved with the data for ease of analysis later
%
%   epochOnsets: trials x onsets vector containing the onset in seconds for
%       the events of interest, relative to time 0 in the trial (from
%       EEG.times)
%
%   epochOffsets: trials x onsets vector containing the offset in seconds
%       for the events of interest, relative to time 0 in the trial (from 
%       EEG.times) 
%
%   saveStem: string that includes the path to the folder to store the file
%       as well as the initial segment of the fileName (the script will
%       append the electrode name); for example
%       '/path/to/power/UCDMC15_TeleporterA_epoched_power_'
%
%   frequencies: vector of log-spaced frequencies to use for analysis
%
%   width: number of cycles to use for wavelet convolution; default = 6
%
%   bufferCycles: number of cycles to use as mirrored buffer during power
%       estimation; default = 3
%

% set wavelet width if not specified
if nargin < 6
    bufferCycles = 3;
    if nargin < 5
        width = 6;
    end
end

% check that we have the same number of onsets and offsets
if (size(epochOnsets) ~= size(epochOffsets))
    error('Different number of onsets and offsets.')
end

% Calculate amount of buffer to use. Based on Mike Cohen's recommendation,
% (p.77 of Analyizng Neural Time Series Data), we'll use the duration of 3
% cycles at the lowest frequency.
bufferMS   = bufferCycles * 1/min(frequencies) * 1000;
bufferBins = bufferMS * 1/1000 * EEG.srate;

% Initialize output array to save
allPowerData = cell(size(epochOnsets, 1), size(epochOnsets, 2));


% Loop over channels
for thisChan = 1:length(chanNames)
    
    % Let the user know the current status
    fprintf(['\n\n\nCalculating power for channel ' num2str(thisChan) ' of ' num2str(length(chanNames)) '\n']);
      
    % find the index of this channel
    chanInd = find(strcmpi(chanNames{thisChan}, {EEG.chanlocs.labels}));
    
    for thisTrial = 1:size(EEG.data, 3)
        
        % Let the user know the current status
        fprintf(['Trial # (' num2str(size(EEG.data,3)) ' total): ' num2str(thisTrial) '\n']);
        
        % Extract data for this trial
        trialData = EEG.data(chanInd, :, thisTrial);
        
        for thisEpoch = 1:size(epochOnsets, 2)
            
            % Find the EEG data within our epoch
            epochInds = find(EEG.times >= epochOnsets(thisTrial, thisEpoch) & EEG.times <= epochOffsets(thisTrial, thisEpoch));
            epochData = trialData(epochInds);
            
            % Create the buffer for the epoched data so that we have enough
            % data to estimate low-frequency power
            numBuffers = ceil(bufferBins / length(epochData));
            bufferData = [];
            
            for thisChunk = 1:numBuffers
                
                % Mirror the signal on every other chunk so that it's
                % continuous
                if mod(thisChunk, 2) == 1
                    bufferData = cat(2, bufferData, fliplr(epochData));
                else
                    bufferData = cat(2, bufferData, epochData);
                end
                
                
            end % thisChunk
            
            % Trim excess data from the buffer
            bufferData = bufferData(1:bufferBins);
            
            % Add the buffer to the epoched data
            epochDataBuffered = [fliplr(bufferData), epochData, bufferData];
            
            % Extract power
            powerVal = single(multienergyvec(epochDataBuffered, frequencies, EEG.srate, width));
            
            % Trim the data to remove the buffers
            powerValTrimmed = powerVal(:, bufferBins + 1 : length(powerVal) - bufferBins);
            
            % Take the log(power) and then z-score it within each frequency
            logPower = log(powerValTrimmed);
            zscorePower = zscore(logPower, 0, 2);
            
            % Add to the summary array
            allPowerData{thisTrial, thisEpoch} = powerValTrimmed;
            
            
        end % thisEpoch
        
    end % thisTrial
    
    % Save the output
    save([saveStem chanNames{thisChan} '.mat'], 'allPowerData', 'trialTypes', 'frequencies');
    
end % thisChan



function [B,t,f]=multienergyvec(S,f,Fs,width)
% function [B,t,f]=multienergyvec(S,f,Fs,width)
% compute the wavelet transform of the input data at the
% frequencies specified by f
% INPUT Args:
% s : signal
% f: vector with frequencies
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).
% OUTPUT Args:
% B: matrix of length(freqs) by length(S) with the wavelet
% coefficients
% t: time vector corresponding to the signal
% f: vector with frequencies for which power was computed
B = zeros(length(f),length(S));

% fprintf('frequency ');
for a=1:length(f)
%     fprintf('%d ',a);
    B(a,:)=energyvec(f(a),S,Fs,width);
end
% fprintf('\n');

t = (1:size(S,2))/Fs;
