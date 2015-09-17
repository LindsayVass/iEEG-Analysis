function signifIntervals = identifySignifIntervals(binaryVec, numContigVal)

% identifySignifIntervals: determine whether there is a string of 1s in the
% binaryVec that is at least numContigVal long
% >> signifIntervals = identifySignifIntervals(binaryVec, numContigVal)
%
% Inputs:
%   binaryVec: vector of binary values (1,0)
%   numContigVal: number of contiguous 1s that must be observed for a
%       significant effect
%
% Output:
%   signifIntervals: 2-column matrix containing the onset (col1) and offset
%       (col2) of each significant interval
%
% Lindsay Vass
% 16 September 2015

% Find the edges where powerBin changes from 0 --> 1 or 1 --> 0
diffBinaryVec = diff(binaryVec);
startInt = find(diffBinaryVec == 1) + 1;
stopInt  = find(diffBinaryVec == -1);

% Handle special edge cases
if (isempty(startInt) && isempty(stopInt))
    if(find(binaryVec) > 0)
        % All values exceed powerThresh
        signifIntervals = [1 length(binaryVec)]; 
    else
        % No values exceed powerThresh
        signifIntervals = [];
    end
elseif (isempty(startInt))
    % Starts on an episode and then stops
    signifIntervals = [1 stopInt];
elseif (isempty(stopInt))
    % Starts an episode and continues until end of EEG
    signifIntervals = [startInt length(binaryVec)];
else
    if (startInt(1) > stopInt(1))
        % We started in an episode
        startInt = [1; startInt];
    end
    if (stopInt(length(stopInt)) < startInt(length(startInt)))
        % We ended on an episode
        stopInt = [stopInt; length(binaryVec)];
    end
    signifIntervals = [startInt stopInt];
end % special edge cases

if (~isempty(signifIntervals)) 
    % Find epochs that exceed the durationBinThresh
    goodEpochs = find((signifIntervals(:, 2) - signifIntervals(:, 1) + 1) >= numContigVal);
    
    if (isempty(goodEpochs))
        signifIntervals = [];
    else
        signifIntervals = signifIntervals(goodEpochs, :);
    end
end