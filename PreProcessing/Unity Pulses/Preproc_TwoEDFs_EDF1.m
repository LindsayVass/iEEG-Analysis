% UCDMC15 Preprocessing for Teleporter A. This dataset exported from UCDMC
% as two EDF files, indicating that there is time missing from the full
% dataset. As a consequence, we'll pre-process each EDF separately and
% combine the data later.
%
% Lindsay Vass 24 April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%% Set up paths
subject_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/';
subject_id  = 'UCDMC15';
edf_file    = 'UCDMC15_04_01_15_teleporter.edf'; % in 'Raw Data' folder
save_stem   = 'UCDMC15_teleporterA_EDF1';
pulsesFile  = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Raw Data/pulses_teleporterA.mat';

addpath(genpath(subject_dir))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/iEEG/Amber Scripts/'))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b/'))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/functions/'))

%% Load patient LFP data in .edf format
cd ([subject_dir '/Raw Data/'])
[hdr,data] = edfread(edf_file);

%% Put in EEG struct
eeglab

EEG.setname = subject_id;
EEG.filename = subject_id;
EEG.srate = hdr.samples(1);

[hour,minute] = strtok(hdr.starttime,'.');
minute = minute(2:end);
[minute,sec] = strtok(minute,'.');
sec = sec(2:end);
EEG.start_time = ([num2str(hour),':',num2str(minute),'.',num2str(sec)]);

% In order to determine the reference time series, you should determine
% which channels are bad and leave them out of the average. Average all
% time series together to get the reference time series.
hdr.label = cellfun(@(x) x(4:end-3),hdr.label(1:end),'UniformOutput',false);

% Work on channel locations, inserting dummy information so we can use ica
% and plot channel properties
ctr=1;

for e = 1:size(data,1)
    
    if isempty(strfind(hdr.label{e},'empty')); % only use real channels
        EEG.data(ctr,:) = data(e,:);
        EEG.chanlocs(ctr).labels = hdr.label{e};
        EEG.chanlocs(ctr).theta = 1;
        EEG.chanlocs(ctr).radius = 1;
        EEG.chanlocs(ctr).X = e;
        EEG.chanlocs(ctr).Y = 1;
        EEG.chanlocs(ctr).Z = 1;
        EEG.chanlocs(ctr).type = e;
        %EEG.ref_vec(ctr) = nkdata.ref_vec(e);
        EEG.chanlocs(ctr).lobe = NaN;
        EEG.chanlocs(ctr).gyrus = NaN;
        EEG.chanlocs(ctr).spikes = NaN;
        EEG.chanlocs(ctr).seizure_onsets = NaN;
        ctr=ctr+1;
    end
    
end

EEG.nbchan = size(EEG.data,1);
eeglab redraw;
pop_eegplot(EEG,1,1,1);

% Save current state of EEG
save([save_stem '_raw.mat'],'EEG');

%% Find marker channel and remove bad channels

% Prompt user to input marker channel number
prompt    = '\n\nInput marker channel number: ';
markerChanNumber = input(prompt);
hdr.label(markerChanNumber) = {'marker'};
EEG.chanlocs(markerChanNumber).labels = 'marker';
eeglab redraw % updates changes in the EEG structure in the GUI

% Prompt user for which channels to remove
badChans = input('\n\nWhich channel(s) should be removed? If multiple, enclose in brackets []: ');
EEG.data(badChans,:) = [];
EEG.nbchan = EEG.nbchan - length(badChans);
EEG.chanlocs([badChans]) = [];
pop_eegplot(EEG,1,1,1);

% Save current state of EEG
save([save_stem '_badChansRemoved.mat'],'EEG');

%% Find the first pulse

% Plot the marker channel data and ask the user to identify the time of the
% first pulse
h = figure;
x = [1:1:size(data,2)]/EEG.srate;
plot(x,data(markerChanNumber,:));
title('Marker channel activity');
xlabel('Time (seconds)');
prompt = '\n\nZoom in on the "Marker channel activity" plot until you can see the first pulse. \nInput the time in seconds just before the first pulse: ';
firstPulseSec = input(prompt);
firstPulseBin = firstPulseSec*EEG.srate;
close(h);

% Remove data before first pulse
bufferSec = 5;
bufferBin = bufferSec*EEG.srate;
maxTrim   = firstPulseBin - bufferBin - 1;

data(:,1:maxTrim) = [];
EEG.data(:,1:maxTrim) = [];
EEG.times(1:maxTrim) = [];
EEG.pnts = size(data,2);
eeglab redraw;
pop_eegplot(EEG,1,1,1);

% Save current state of EEG
save([save_stem '_badChansRemoved_trimmed.mat'],'EEG');

%% Find all remaining pulses

% Update time of the first pulse
firstPulseBin = firstPulseBin - maxTrim - 1;
firstPulseSec = bufferSec;

% Plot time around the first pulse
plotBufferBin = 250;
markerData    = data(markerChanNumber,:);
h = figure;
ax = gca;
plot(markerData(firstPulseBin:firstPulseBin+(plotBufferBin*2)));
title('First pulse?');

% Ask user if the current plot includes the first pulse. Keep shifting
% ahead in time until it does.
prompt = 'Does the "First Pulse?" plot contain the first pulse?';
answer = questdlg(prompt);

numTries = 1;
while 1
    if (strcmpi(answer,'Yes') == 1)
        break
    end
    
    firstPulseBin = firstPulseBin + plotBufferBin*2;
    numTries = numTries + 1;
    
    plot(ax, markerData(firstPulseBin-plotBufferBin:firstPulseBin+plotBufferBin));
    title('First pulse?');
    answer = questdlg(prompt);
    
end

% Get pulse polarity (positive/negative)
prompt = 'Is the pulse polarity positive or negative?';
pulsePolarity = questdlg(prompt,'Pulse Polarity','Positive','Negative','Positive');

% Find the max value within the range and plot
if (numTries == 1)
    plotRange = markerData(firstPulseBin:firstPulseBin + (plotBufferBin * 2));
else
    plotRange = markerData(firstPulseBin-plotBufferBin:firstPulseBin+plotBufferBin);
end

if(strcmpi(pulsePolarity,'Positive') == 1)
    maxVal = max(plotRange);
    maxValInd = find(plotRange == maxVal);
    if(length(maxValInd) > 1)
        maxValInd = maxValInd(1);
    end
    y = get(ax,'ylim');
    hold on;
    plot(ax,[maxValInd maxValInd],y,'--r');
else
    minVal = min(plotRange);
    minValInd = find(plotRange == minVal);
    if(length(minValInd) > 1)
        minValInd = minValInd(1);
    end
    y = get(ax,'ylim');
    hold on;
    plot(ax,[minValInd minValInd],y,'--r');
end

% Ask user to check that the peak of the pulse was found correctly
prompt = 'Does the dotted red line mark the peak of the pulse?';
answer = questdlg(prompt);

if(strcmpi(answer,'Yes') == 1)
    if(strcmpi(pulsePolarity,'Positive') == 1)
        
        firstPulseBin = firstPulseBin - plotBufferBin + maxValInd;
        close(h);
        
    else
        firstPulseBin = firstPulseBin - plotBufferBin + minValInd;
        close(h)
    end
else
    disp('Peak of pulse not automatically found. QUITTING!');
    forceError;
end

%% Find all pulses

load(pulsesFile);

% Pulses are ticks since Jan 1, 0001. Re-reference them to the first pulse.
pulsesRef = pulses - pulses(1);

% Convert from ticks to bins
ticksPerSec = 10000000;
ticksPerBin = ticksPerSec/EEG.srate;
pulsesBin   = round(pulsesRef/ticksPerBin) + firstPulseBin;

% Ask for visual confirmation of each pulse
prompt = 'Does the plot show a pulse with a dotted red line at the peak?';
pulsesBinConf = firstPulseBin;
currBin = firstPulseBin + EEG.srate;

pulseNum = 2;

while 1
    
    h = figure;
    ax = gca;
    
    if (currBin > length(markerData))
        break
    end
    
    x = currBin-plotBufferBin:currBin+plotBufferBin;
%     x = currBin-plotBufferBin:currBin+plotBufferBin*3; % use if problem is that pulse is delayed
    y = markerData(x);
    
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
    peakDiff = currBin - x(m);
    
    answer = MFquestdlg([0.5, 0.5],prompt);
    
    if(strcmpi(answer,'Yes') == 1)
        pulsesBinConf(end+1) = x(m);
        pulseNum = pulseNum + 1;
        currBin = currBin + EEG.srate - peakDiff;
        close(h);
    else
        y2 = y';
        fprintf('\n\nThere are 2 ways the peak finder could fail. First, the entire pulse might be \ndelayed. In this case, increase the range for ?x? until the pulse is visible \nin the plot. Second, the pulse might have a strange shape, with multiple \npeaks. In this case, examine the y2 variable and find the local peak that \ncorresponds to the first peak of the pulse. Then set m = peak. For example, \nif your desired peak occurs at index 250 of y2, then m = 250. Then, execute \nthe remainder of the loop (line 245 onwards). If everything is good, execute \nthe while loop again and it will pick up where you left off.\n\n');
        
        return;
    end
    
end

% Save the confirmed pulses
indEEG = pulsesBinConf';
unityTicks = pulses(1:length(indEEG));
save([save_stem '_pulse_timing.mat'],'indEEG','unityTicks');