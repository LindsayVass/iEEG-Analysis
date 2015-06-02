function binaryVector = findPepisodes(powerVal, powerThresh, durationThresh, frequency, samplingRate)

% function binaryVector = findPepisodes(powerVal, powerThresh, durationThresh, frequency, samplingRate)
% Identify for each time point whether a sustained oscillation is present
% and return a binary vector that contains a 1 for timepoints with
% oscillations present and a 0 for timepoints without oscillations present.
% To be counted as an oscillation, the power must exceed the power
% threshold and duration threshold.
%
% INPUTS:
%   powerVal: vector of power values to analyze
%   powerThresh: threshold that power must exceed to count as an
%       oscillation
%   durationThresh: number of cycles for which the power must exceed the
%       powerThresh to be counted as an oscillation
%   frequency: frequency in Hz at which pepisode is being detected
%   samplingRate: sampling rate in Hz of the input EEG data
%
% OUTPUT:
%   binaryVector: binary vector that indicates for each timepoint whether a
%       sustained oscillation was present (1) or not (0)

if nargin < 5
    error ('Required arguments not supplied.')
end

% Initialize output vector
binaryVector = zeros(1, length(powerVal));

% Convert durationThresh from cycles to bins
durationBinThresh = durationThresh * samplingRate / frequency;

% Zero all values under the powerThresh
powerVal(find(powerVal<powerThresh)) = 0;

% Binary version of thresholded power data
powerBin = (powerVal > 0);

% Find the edges where powerBin changes from 0 --> 1 or 1 --> 0
diffPowerBin = diff(powerBin);
startOsc = find(diffPowerBin == 1) + 1;
stopOsc  = find(diffPowerBin == -1) + 1;

% Handle special edge cases
if (isempty(startOsc) && isempty(stopOsc))
    if(find(powerBin) > 0)
        % All values exceed powerThresh
        aboveThreshInterval = [1; length(powerVal)]; 
    else
        % No values exceed powerThresh
        aboveThreshInterval = [];
    end
elseif (isempty(startOsc))
    % Starts on an episode and then stops
    aboveThreshInterval = [1; stopOsc];
elseif (isempty(stopOsc))
    % Starts an episode and continues until end of EEG
    aboveThreshInterval = [startOsc; length(powerVal)];
else
    if (startOsc(1) > stopOsc(1))
        % We started in an episode
        startOsc = [1 startOsc];
    end
    if (stopOsc(length(stopOsc)) < startOsc(length(startOsc)))
        % We ended on an episode
        stopOsc = [stopOsc length(powerVal)];
    end
    aboveThreshInterval = [startOsc; stopOsc];
end % special edge cases

if (~isempty(aboveThreshInterval)) 
    % Find epochs that exceed the durationBinThresh
    goodEpochs = find((aboveThreshInterval(2,:) - aboveThreshInterval(1,:)) >= durationBinThresh);
    
    if (isempty(goodEpochs))
        aboveThreshInterval = [];
    else
        aboveThreshInterval = aboveThreshInterval(:, goodEpochs);
    end
    
    % Make the binary vector for output
    for thisGoodEpoch = 1:size(aboveThreshInterval, 2)
        binaryVector(aboveThreshInterval(1, thisGoodEpoch):aboveThreshInterval(2, thisGoodEpoch)) = 1;
    end
end


end