function [onsetTime, offsetTime, duration] = calcPepisodeSustainedness(powerVal, powerThresh, durationThresh, frequency, samplingRate, eegTimes)

% USAGE: onsetTime, offsetTime, duration] = calcPepisodeSustainedness(powerVal, powerThresh, durationThresh, frequency, samplingRate)
%
% Identify sustained oscillations using the supplied amplitude and duration
% thresholds. Then, for each oscillatory episode, return the timing of the
% onset and offset as well as the duration.
%
% INPUTS:
%   powerVal: vector of power values to analyze
%   powerThresh: threshold that power must exceed to count as an
%       oscillation
%   durationThresh: number of cycles for which the power must exceed the
%       powerThresh to be counted as an oscillation
%   frequency: frequency in Hz at which pepisode is being detected
%   samplingRate: sampling rate in Hz of the input EEG data
%   eegTimes: times (in seconds) of each bin of the epoch (can be found in
%       EEG.times of input data structure)
%
% OUTPUT:
%   onsetTime: time at which the oscillation started, where time 0 is
%       the time of teleporter entry
%   offsetTime: time at which the oscillation ended, where time 0 is
%       the time of teleporter entry
%   duration: length of time that the oscillation was sustained
%       (i.e., offsetTime - onsetTime)


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

if (size(eegTimes, 1) == 1)
    eegTimes = eegTimes';
end

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
        onsetTime           = [];
        offsetTime          = [];
        duration            = [];
    else
        aboveThreshInterval = eegTimes(aboveThreshInterval(:, goodEpochs));
        % Make the output
        for thisGoodEpoch = 1:size(aboveThreshInterval, 2)
            onsetTime(thisGoodEpoch)  = aboveThreshInterval(1, thisGoodEpoch);
            offsetTime(thisGoodEpoch) = aboveThreshInterval(2, thisGoodEpoch);
            duration(thisGoodEpoch)   = aboveThreshInterval(2, thisGoodEpoch) - aboveThreshInterval(1, thisGoodEpoch);
        end
    end
else
    onsetTime           = [];
    offsetTime          = [];
    duration            = [];
end


end