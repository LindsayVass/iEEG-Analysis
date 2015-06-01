function powerDistribution = calcPowerDistributionLKV(eegList, chanName, frequencies, waveletWidth, saveFile)

% calcPowerDistrubtionLKV Calculate distribution of power at a given set of
% frequencies using Morlet wavelets.
%
% powerDistribution = calcPowerDistributionLKV(eegList, chanNames,
% frequencies, waveletWidth, saveStem)
%
% INPUTS:
%   eegList: cell array of paths to the EEGLAB dataset(s) on which to
%       calculate power
%   chanName: string of electrode name to use for analysis
%   frequencies: vector of frequencies in Hz at which to calculate power
%       distributions
%   waveletWidth (optional): number of Morlet wavelet cycles to use for
%       power estimation; default = 6
%   saveFile (optional): string representing the path of the desired
%       power distribution output file; example =
%       '/path/to/distrib/UCDMC13_TeleporterA_LAD1_power_distribution.mat'
%
% OUTPUT:
%   powerDistribution: structure with fields for the frequency at which the
%       distribution was calculated and the distribution itself
%
% Lindsay Vass
% v1: 1 June 2015

% Specify waveletWidth if not indicated by user
if nargin < 4
    warning('Argument waveletWidth not specified. Using default value 6.');
    waveletWidth = 6;
end

if nargin < 3
    error('Arguments eegList, chanNames, and frequencies are required.')
end

% Make sure the EEG files exist
for thisFile = 1:length(eegList)
    if ~exist(eegList{thisFile}, 'file')
        error(['EEG file does not exist: ' eegList{thisFile}])
    end
end

% Initialize output struct
powerDistribution = struct;

% Initialize struct to hold all of the power values from each segment of
% each EEG
powerHolder = struct;

% Loop through EEG datasets
for thisEEG = 1:length(eegList)
    
    % Load the EEG data
    EEG = pop_loadset(eegList{thisEEG});
    
    % Identify the boundaries within the EEG file
    boundaries = findBoundaries(EEG);
    
    % find the index of this channel
    chanInd = find(strcmpi(chanName, {EEG.chanlocs.labels}));
    
    % Calculate power
    [eegPower] = extractPowerWithBoundaries(EEG, chanInd, boundaries, frequencies, waveletWidth);
    
    % Add to summary array
    powerHolder(thisEEG).EEG = eegPower;
    
end % thisEEG

% Build the power distribution by concatenating all power values
allPower = [];

for thisEEG = 1:size(powerHolder, 2)
    
    for thisSegment = 1:size(powerHolder(thisEEG).EEG, 2)
        
        allPower = cat(2, allPower, powerHolder(thisEEG).EEG(thisSegment).segment.power);
        
    end % thisSegment
    
end % thisEEG

% If there's zero power, this will cause an error when we take the log, so
% remove them
allPower(allPower == 0) = NaN;

% Take the log of the power values
logPowerVal = log10(double(allPower));

% Get the mean log(Power) for each frequency
logPowerMean = nanmean(logPowerVal, 2);

% Fit a chi-square distribution to the power values across
% frequencies and return the cumulative distribution function
% (cdf). The first value of the cdf are the frequencies, and
% the remaining 1000 values are the distribution.
[cdf, ~] = chi_squarefit(frequencies, logPowerMean);
cdf = cdf(2:end, :)';

% Loop through frequencies and add data to output struct
for thisFreq = 1:length(frequencies)
    
    powerDistribution(thisFreq).frequency = frequencies(thisFreq);
    powerDistribution(thisFreq).distribution = cdf(thisFreq, :);
    
end % thisFreq

if exist('saveFile', 'var')
    
    save(saveFile, 'powerDistribution');
    
end


end



function [thresh,R2] = chi_squarefit(freqs, pows)
% function [thresh,R2] = chi_squarefit(freqs, pows)
% calculate the fit for the pepisode
nbins=1000;

pv=polyfit(log10(freqs),pows',1); % linear regression
means=10.^(polyval(pv,log10(freqs))); % these are the fitted mean powers

R=corrcoef(polyval(pv,log10(freqs)),pows');
R2=R(1,2)^2;

% Now, compute the chi2cdf vals.
thresh=[freqs;chi2inv(0:(1/nbins):(1-(1/nbins)),2)'*(means/2)];
% divide the means by two because the mean of the chi-square distribution is equal to the computed mean divided by the degrees of freedom

end


