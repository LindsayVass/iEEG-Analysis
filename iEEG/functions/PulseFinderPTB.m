function [pulsesBinConf] = PulseFinderPTB(EEG, markerChanNumber, startBin, direction, eventList, plotBuffer)
% Find pulses in EEG data. Will save to the subject's "Mat Files"
% directory.
%
% [indEEG, unityTicks] = PulseFinder(subject_dir, subject_id, teleporter,
% EEG, markerChanNumber, startBin, direction, plotBuffer)
%
% INPUTS:
%   subject_dir: main directory for subject
%   subject_id: ID for subject (e.g., UCDMC15)
%   teleporter: which teleporter (e.g., TeleporterA)
%   EEG: EEG structure to be used for pulse finding
%   markerChanNumber: which channel of EEG is the marker
%   startBin: index of the first or last pulse of the session
%   direction: 'forward' or 'backward' depending on if it's EDF1 or EDF2
%   eventList: structure output by psychtoolbox script
%   plotBuffer: how many bins on either side of the pulse to plot; default = 250;
%
% OUTPUT:
%   pulsesBinConf: vector of the peaks of the pulses in EEG bins


if nargin == 5
    plotBuffer = 250;
end

% initialize vector to hold confirmed pulses
pulsesBinConf = [];

pulseNum = 1;

% set direction of pulse finding
if strcmpi(direction,'forward') == 1
    pulseDir = 1;
elseif strcmpi(direction, 'backward') == 1
    pulseDir = -1;
else
    error('Not a valid direction. Use either "forward" or "backward".');
end


markerData = EEG.data(markerChanNumber,:);
h = figure;
ax = gca;

%% Plot the starting pulse
if pulseDir == 1
    
    x = startBin:startBin+(plotBuffer*2);
    y = markerData(x);
    
    plot(ax, x, y);
    title('First pulse?');
    
    % Ask user if the current plot includes the first pulse. Keep shifting
    % ahead in time until it does.
    prompt = 'Does the "First Pulse?" plot contain the first pulse?';
    answer = questdlg(prompt);
    
    while 1
        if (strcmpi(answer,'Yes') == 1)
            break
        end
        
        startBin = startBin + plotBuffer*2;
        x = startBin:startBin+(plotBuffer*2);
        y = markerData(x);
        
        plot(ax, x, y);
        title('First pulse?');
        answer = questdlg(prompt);
        
    end
else
    x = startBin-plotBuffer:startBin+plotBuffer;
    y = markerData(x);
    
    plot(ax, x, y);
    title('Last pulse?');
    
    % Ask user if the current plot includes the last pulse. Keep shifting
    % ahead in time until it does.
    prompt = 'Does the "Last Pulse?" plot contain the last pulse?';
    answer = questdlg(prompt);
    
    while 1
        if (strcmpi(answer,'Yes') == 1)
            break
        end
        
        startBin = startBin - plotBuffer*2;
        x = startBin-plotBuffer:startBin+plotBuffer;
        y = markerData(x);
        
        plot(ax, x, y);
        title('Last pulse?');
        answer = questdlg(prompt);
        
    end
end

%% Get pulse polarity (positive/negative)
prompt = 'Is the pulse polarity positive or negative?';
pulsePolarity = questdlg(prompt,'Pulse Polarity','Positive','Negative','Positive');

% Find the max value within the range and plot
% x = startBin:startBin+plotBuffer*2;
% y = markerData(x);

if(strcmpi(pulsePolarity,'Positive') == 1)
    maxVal = max(y);
    maxValInd = find(y == maxVal);
    if(length(maxValInd) > 1)
        maxValInd = maxValInd(1);
    end
    ylims = get(ax,'ylim');
    hold on;
    plot(ax,[x(maxValInd) x(maxValInd)],ylims,'--r');
else
    minVal = min(y);
    minValInd = find(y == minVal);
    if(length(minValInd) > 1)
        minValInd = minValInd(1);
    end
    ylims = get(ax,'ylim');
    hold on;
    plot(ax,[x(minValInd) x(minValInd)],ylims,'--r');
end

%% Ask user to check that the peak of the pulse was found correctly
prompt = 'Does the dotted red line mark the peak of the pulse?';
answer = questdlg(prompt);

if(strcmpi(answer,'Yes') == 1)
    if(strcmpi(pulsePolarity,'Positive') == 1)
        
        startBin = x(maxValInd);
        close(h);
        
    else
        startBin = x(minValInd);
        close(h)
    end
else
    expectedBin = startBin;
    FindOtherPeaks;
end

%% Find all remaining pulses

% Ask for visual confirmation of each pulse
prompt = 'Does the plot show a pulse with a dotted red line at the peak?';
pulsesBinConf = nan(length(eventList.time),1);
pulsesBinConf(1) = startBin;

peakDiff = 0;
previousPulseBin = startBin;
for thisPulse = 2:length(pulsesBinConf)
    
    % calculate the interval in seconds between this pulse and the previous
    % one
    pulseInterval = etime(eventList.time{thisPulse}, eventList.time{thisPulse - 1});
    
    % convert this interval to EEG bin
    expectedBin = round(pulseInterval * EEG.srate) + previousPulseBin - peakDiff;
%     expectedBin = round(pulseInterval * EEG.srate) + previousPulseBin;
    
    % plot the marker data around the expected time of the pulse
    h = figure;
    ax = gca;
    x = expectedBin - plotBuffer:expectedBin + plotBuffer;
    y = markerData(x);
    
    % Make sure this pulse is within the data range
    if pulseDir == 1
        if (expectedBin + plotBuffer > length(markerData))
            break
        end
    else
        if (expectedBin - plotBuffer < 1)
            break
        end
    end
    
    % find local peak
    if(strcmpi(pulsePolarity,'Positive') == 1)
        m = find(y == max(y));
    else
        m = find(y == min(y));
    end
    
    % sometimes 2 adjacent time points may have the same value. If this
    % happens, take the first.
    if length(m) > 1
        m = m(1);
    end
    
    plot(ax, x, y);
    yLims = get(ax,'ylim');
    hold on;
    plot(ax, [x(m) x(m)], yLims, '--r');
    hold off
    
    % get difference between the expected peak (currBin) and the actual
    % peak so we can account for drift
    peakDiff = expectedBin - x(m);
    
    answer = MFquestdlg([0.5, 0.5],prompt);
    
    if(strcmpi(answer,'Yes') == 1) % go to the next pulse
        pulsesBinConf(thisPulse) = x(m);
        previousPulseBin = x(m);
        close(h);
        
    else % something's wrong with the pulse
        
        question = 'What is the problem with the pulse?';
        questTitle = 'Pulse Debug';
        button1 = 'There is not a pulse';
        button2 = 'The peak is not marked correctly';
        answer = questdlg(question, questTitle, button1, button2, button2);
        
        % If the pulse was delayed, move forward until we
        % find it
        if strcmpi(answer, button1) == 1
            prompt = 'Is there a pulse shown now?';
            while 1
                
                if pulseDir == 1
                    expectedBin = expectedBin + plotBuffer*2;
                else
                    expectedBin = expectedBin - plotBuffer*2;
                end
                x = expectedBin-plotBuffer:expectedBin+plotBuffer;
                y = markerData(x);
                
                plot(ax, x, y);
                title('Is there a pulse');
                answer = questdlg(prompt);
                
                if (strcmpi(answer,'Yes') == 1)
                    
                    % find local peak
                    if(strcmpi(pulsePolarity,'Positive') == 1)
                        m = find(y == max(y));
                    else
                        m = find(y == min(y));
                    end
                    
                    % sometimes 2 adjacent time points may have the same value. If this
                    % happens, take the first.
                    if length(m) > 1
                        m = m(1);
                    end
                    
                    plot(ax, x, y);
                    yLims = get(ax,'ylim');
                    hold on;
                    plot(ax, [x(m) x(m)], yLims, '--r');
                    hold off
                    
                    % get difference between the expected peak (currBin) and the actual
                    % peak so we can account for drift
                    peakDiff = expectedBin - x(m);
                    prompt = 'Does the plot show a pulse with a dotted red line at the peak?';
                    answer = MFquestdlg([0.5, 0.5],prompt);
                    
                    if(strcmpi(answer,'Yes') == 1)
                        pulsesBinConf(thisPulse) = x(m);
                        previousPulseBin = x(m);
                        close(h);
                        break
                    else
                        m = FindOtherPeaks(pulsePolarity, ax, x, y);
                        pulsesBinConf(thisPulse) = x(m);
                        previousPulseBin = x(m);
                        close(h);
                        break
                    end
                end
            end
        else
            m = FindOtherPeaks(pulsePolarity, ax, x, y);
            peakDiff = expectedBin - x(m);
            pulsesBinConf(thisPulse) = x(m);
            previousPulseBin = x(m);
            close(h);
        end
    end
end

end

function [m] = FindOtherPeaks(pulsePolarity,ax,x,y)
% If we selected the wrong peak, find all local
% maxima and ask for the correct one
prompt = 'Is the peak properly marked now?';

if strcmpi(pulsePolarity,'Positive') == 1
    [~, peakList] = findpeaks(double(y));
else
    yInv = 1.01*max(y) - y;
    [~, peakList] = findpeaks(double(yInv));
end

% loop through peaks and highlight each one in
% turn
thisPeak = 1;
while 1
    
    plot(ax, x, y);
    hold on;
    ylims = get(gca, 'YLim');
    plot([x(peakList(thisPeak)) x(peakList(thisPeak))],ylims, '--r');
    hold off;
    
    title('Is the peak properly marked?');
    answer = questdlg(prompt);
    
    if (strcmpi(answer,'Yes') == 1)
        m = peakList(thisPeak);
        prompt = 'Does the plot show a pulse with a dotted red line at the peak?';
        break
    else
        thisPeak = thisPeak + 1;
        if thisPeak > length(peakList)
            error('Out of peaks to try. Quitting now.')
        end
    end
    
end
end