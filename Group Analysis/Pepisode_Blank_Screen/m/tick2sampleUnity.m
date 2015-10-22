function [thisSample, EDF] = tick2sampleUnity(thisTick, EDF1ticks, EDF2ticks, EDF1samples, EDF2samples)
% Convert from ticks to EEG samples using Unity pulses
% >> thisSample = tick2sampleUnity(thisTick, EDF1ticks, EDF2ticks, EDF1samples, EDF2samples)
%
% Inputs:
%   thisTick: time (in ticks) to convert from
%   EDF1ticks: list of tick times for EDF1
%   EDF2ticks: list of tick times for EDF2
%   EDF1samples: list of sample times for EDF1
%   EDF2samples: list of sample times for EDF2
%
% Output:
%   thisSample: time (in samples) of thisTick
%   EDF: which EDF is thisTick found in (1 or 2)

if thisTick - EDF1ticks(end) < 0 % if in EDF1
    
    thisSample = fitTicks(thisTick, EDF1ticks, EDF1samples, 5);
    EDF = 1;
    
elseif thisTick - EDF2ticks(1) <0 % trial in lost EEG between EDFs
    
    thisSample = NaN;
    EDF = NaN;
    
else % in EDF2
    
    thisSample = fitTicks(thisTick, EDF2ticks, EDF2samples, 5);
    EDF = 2;
    
end


function thisSample = fitTicks(thisTick, tickList, sampleList, fitRange)

% find the pulse time in ticks closest to thisTick
tickDiff   = abs(thisTick - tickList);
minDiffInd = find(tickDiff == min(tickDiff));

% get a range of values before and after
fitInds = [minDiffInd - fitRange:1:minDiffInd + fitRange];
fitInds(fitInds < 1) = [];
fitInds(fitInds > length(tickList)) = [];

% get fit of line using the inds selected above
fitP = polyfit(tickList(fitInds), sampleList(fitInds), 1);
thisSample = round(thisTick * fitP(1) + fitP(2));