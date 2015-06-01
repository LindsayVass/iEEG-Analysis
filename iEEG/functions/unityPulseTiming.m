function [ticks EEG1bins EEG2bins] = unityPulseTiming(eventTimingInTicks, EDF1PulseFile, EDF2PulseFile, fitRange)
% Converts timing in ticks to EEG bins using the provided pulse files. This
% function will identify the pulses surrounding the time point of interest
% and fit a linear model to them to convert from timing in ticks to timing
% in EEG indices. This function will display the linear fit to the user for
% confirmation.
%
% OUTPUTS:
% ticks - the vector of ticks provided by eventTimingInTicks
% EEGbins - vector of EEG indices for values in ticks
%
% Required arguments:
% eventTimingInTicks - vector of tick values
% EDF1PulseFile - mat file containing two variables, indEEG and unityTicks,
% that contain the mapping between pulses as ticks and pulses as EEG
% indices
%
% Optional arguments:
% EDF2PulseFile - optional second pulse timing file
% fitRange - # of ticks before and after time point to use for linear
% model; default = 5


if ~exist('eventTimingInTicks', 'var') || ~exist('EDF1PulseFile', 'var')
    error('eventTimingInTicks and EDF1PulseFile are required inputs');
end

if ~exist('fitRange', 'var')
    fitRange = 5;
end

% initialize output arrays
ticks = eventTimingInTicks;
EEG1bins = [];
EEG2bins = [];

% load pulse file(s)
load(EDF1PulseFile);
unityTicks1 = unityTicks;
indEEG1     = indEEG;

if exist('EDF2PulseFile', 'var')
    load(EDF2PulseFile);
    unityTicks2 = unityTicks;
    indEEG2     = indEEG;
end

for thisEvent = 1:length(eventTimingInTicks)
    
    % select the current tick
    thisTick = eventTimingInTicks(thisEvent);
    
    if (thisTick - unityTicks1(end) < 0) % if trial in EDF1
        
        % find the pulse time (in ticks) closest to our tick
        tickDiff  = abs(thisTick - unityTicks1);
        minDiffInd = find(tickDiff == min(tickDiff));
        
        % get a range of values before and after
        fitInds = [minDiffInd-fitRange:1:minDiffInd+fitRange];
        
        % get fit of line using the inds selected above
        fit_P = polyfit(unityTicks1(fitInds),indEEG1(fitInds),1);
        fit_y = polyval(fit_P,unityTicks1(fitInds));
        eventBin = round(thisTick*fit_P(1) + fit_P(2));
        
        h = figure;
        plot(unityTicks1(fitInds),indEEG1(fitInds),'k*')
        hold on;
        plot(unityTicks1(fitInds),fit_y,'-')
        scatter(thisTick,eventBin,'ro')
        
        % make sure the fit looks good or else quit
        answer = questdlg('Does the fit look good?');
        if(strcmpi(answer,'Yes') ~= 1)
            break
        end
        
        close(h);
        
        EEG1bins(end+1) = eventBin;
        
    elseif (thisTick - unityTicks2(1) < 0) % trial in lost EEG between EDF files
        warning('Event not found in time range of EDF file(s) provided');
        continue
    else % trial in EDF2
        
        % find the pulse time (in ticks) closest to our tick
        tickDiff  = abs(thisTick - unityTicks2);
        minDiffInd = find(tickDiff == min(tickDiff));
        
        % get a range of values before and after
        fitInds = [minDiffInd-fitRange:1:minDiffInd+fitRange];
        
        % get fit of line using the inds selected above
        fit_P = polyfit(unityTicks2(fitInds),indEEG2(fitInds),1);
        fit_y = polyval(fit_P,unityTicks2(fitInds));
        eventBin = round(thisTick*fit_P(1) + fit_P(2));
        
        h = figure;
        plot(unityTicks2(fitInds),indEEG2(fitInds),'k*')
        hold on;
        plot(unityTicks2(fitInds),fit_y,'-')
        scatter(thisTick,eventBin,'ro')
        
        % make sure the fit looks good or else quit
        answer = questdlg('Does the fit look good?');
        if(strcmpi(answer,'Yes') ~= 1)
            break
        end
        
        close(h);
        
        EEG2bins(end+1) = eventBin;
        
    end
end

end