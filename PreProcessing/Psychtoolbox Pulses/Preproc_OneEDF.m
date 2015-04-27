% UCDMC14 Preprocessing for Teleporter B. 
% Lindsay Vass 27 April 2015

clear all; close all; clc;

%% Set up paths
subject_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/';
subject_id  = 'UCDMC14';
edf_file    = 'UCDMC14_020715.edf'; % in 'Raw Data' folder
save_stem   = 'UCDMC14_TeleporterB_';

addpath(genpath(subject_dir))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/iEEG/Amber Scripts/'))
addpath(genpath('/Users/Lindsay/Documents/MATLAB/eeglab13_3_2b/'))

%% UDIs
epoch_savefile = ['epoch_teleporterB_' subject_id];

%% Choose behavioral data file
load '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Raw Data/UCDMC14_020715_pre_eyemoves/SubjectUCDMC14_eyemoves_Data_interim.mat';
pretest_master_event_list = master_event_list;

load '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Raw Data/UCDMC14_020715_post/SubjectUCDMC14_2715_post_Data_interim.mat/';
posttest_master_event_list = master_event_list;

%% Load patient LFP data in .edf format
cd ([subject_dir '/Raw Data/'])
[hdr,data] = edfread(edf_file);

%% Put in EEG struct
eeglab

EEG.setname = 'UCDMC14';
EEG.filename = 'UCDMC14';
EEG.srate = hdr.samples(1);

[hour,minute] = strtok(hdr.starttime,'.');
minute = minute(2:end);
[minute,sec] = strtok(minute,'.');
sec = sec(2:end);
EEG.start_time = ([num2str(hour),':',num2str(minute),'.',num2str(sec)]);

% In order to determine the reference time series, you should determine
% which channels are bad and leave them out of the average. Average all
% time series together to get the reference time series.
%EEG.ref_tseries = nkdata.ref_tseries;
hdr.label = cellfun(@(x) x(4:end-3),hdr.label(1:end),'UniformOutput',false);
% hdr.label{hdr.ns} = 'marker';

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
eeglab redraw % updates changes in the EEG structure in the GUI

pop_eegplot(EEG,1,1,1);
save('shortcut','EEG');

%% Find marker channel and remove bad channels

% Prompt user to input marker channel number
prompt    = '\n\nInput marker channel number: ';
markerChanNumber = input(prompt);
hdr.label(markerChanNumber) = {'marker'};
EEG.chanlocs(markerChanNumber).labels = 'marker';
eeglab redraw % updates changes in the EEG structure in the GUI
pop_eegplot(EEG,1,1,1);

% Prompt user for which channels to remove
badChans = input('\n\nWhich channel(s) should be removed? If multiple, enclose in brackets []: ');
EEG.data(badChans,:) = [];
EEG.nbchan = EEG.nbchan - length(badChans);
EEG.chanlocs([badChans]) = [];
pop_eegplot(EEG,1,1,1);

% Save current state of EEG
save([save_stem '_badChansRemoved.mat'],'EEG');


%% FIND THE TIME OF THE FIRST PULSE FOR BOTH PRE-TEST AND POST-TEST
keyboard;

% PRE-TEST: 135
% POST-TEST: 2373

%% Estimate time of first pulse and trim data before it

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


%% Find and characterize polarity of first pulse

% Update time of the first pulse
firstPulseBin = firstPulseBin - maxTrim - 1;
firstPulseSec = bufferSec;

% Plot time around the first pulse
plotBufferBin = 256;
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

%% Get pulse timing from PsychToolbox output

firstPulsePTB = pretest_master_event_list.time{1};

% Initialize vector to hold the pulse timing information and add the first
% pulse to it
prePulseBins = nan(length(pretest_master_event_list.time),1);
prePulseBins(1) = firstPulseBin;

% increase plotBufferBin
plotBufferBin = 512;

prompt = 'Does the plot show a pulse with a dotted red line at the peak?';

previousPulseBin = firstPulseBin;
for thisPulse = 2:length(prePulseBins) % loop through pulses
    
    % calculate the interval in seconds between this pulse and the previous
    % one
    pulseInterval = etime(pretest_master_event_list.time{thisPulse}, pretest_master_event_list.time{thisPulse - 1});
    
    % convert this interval to EEG bin
    expectedBin = round(pulseInterval * EEG.srate) + previousPulseBin;
    
    % plot the marker data around the expected time of the pulse
    h = figure;
    ax = gca;
    x = expectedBin - plotBufferBin:expectedBin + plotBufferBin;
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
    
    answer = MFquestdlg([0.5, 0.5],prompt);
    
    if(strcmpi(answer,'Yes') == 1)
        prePulseBins(thisPulse) = x(m);
        previousPulseBin = x(m);
        close(h);
    else
        y2 = y';
        fprintf('\n\nThere are 2 ways the peak finder could fail. First, the entire pulse might be \ndelayed. In this case, increase the range for x until the pulse is visible \nin the plot. Second, the pulse might have a strange shape, with multiple \npeaks. In this case, examine the y2 variable and find the local peak that \ncorresponds to the first peak of the pulse. Then set m = peak. For example, \nif your desired peak occurs at index 250 of y2, then m = 250. Then, execute \nthe remainder of the loop. If everything is good, update \nthe thisPulse variable and continue the rest of the loop.\n\n');
        
        return;
    end
    
end

%% Estimate time of first pulse of the post-test pulses

% Plot the marker channel data and ask the user to identify the time of the
% first pulse
h = figure;
x = [1:1:size(data,2)]/EEG.srate;
plot(x,data(markerChanNumber,:));
title('Marker channel activity');
xlabel('Time (seconds)');
prompt = '\n\nZoom in on the "Marker channel activity" plot until you can see the first pulse of the post-test pulses. \nInput the time in seconds just before the first pulse: ';
firstPulseSec = input(prompt);
firstPulseBin = firstPulseSec*EEG.srate;
close(h);

%% Get precise timing of first pulse

% Plot time around the first pulse
plotBufferBin = 256;
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

%% Get pulse timing from PsychToolbox output

firstPulsePTB = posttest_master_event_list.time{1};

% Initialize vector to hold the pulse timing information and add the first
% pulse to it
postPulseBins = nan(length(posttest_master_event_list.time),1);
postPulseBins(1) = firstPulseBin;

% increase plotBufferBin
plotBufferBin = 512;

prompt = 'Does the plot show a pulse with a dotted red line at the peak?';

previousPulseBin = firstPulseBin;
for thisPulse = 2:length(postPulseBins) % loop through pulses
    
    % calculate the interval in seconds between this pulse and the previous
    % one
    pulseInterval = etime(posttest_master_event_list.time{thisPulse}, posttest_master_event_list.time{thisPulse - 1});
    
    % convert this interval to EEG bin
    expectedBin = round(pulseInterval * EEG.srate) + previousPulseBin;
    
    % plot the marker data around the expected time of the pulse
    h = figure;
    ax = gca;
    x = expectedBin - plotBufferBin:expectedBin + plotBufferBin;
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
    
    answer = MFquestdlg([0.5, 0.5],prompt);
    
    if(strcmpi(answer,'Yes') == 1)
        postPulseBins(thisPulse) = x(m);
        previousPulseBin = x(m);
        close(h);
    else
        y2 = y';
        fprintf('\n\nThere are 2 ways the peak finder could fail. First, the entire pulse might be \ndelayed. In this case, increase the range for x until the pulse is visible \nin the plot. Second, the pulse might have a strange shape, with multiple \npeaks. In this case, examine the y2 variable and find the local peak that \ncorresponds to the first peak of the pulse. Then set m = peak. For example, \nif your desired peak occurs at index 250 of y2, then m = 250. Then, execute \nthe remainder of the loop. If everything is good, update \nthe thisPulse variable and continue the rest of the loop.\n\n');
        
        return;
    end
    
end

%% Perform regression for syncing EEG rig to testing laptop
% Get the timestamps from the stimulus computer for both sets of pulses
pre_time_orig = pretest_master_event_list.time;
post_time_orig = posttest_master_event_list.time;

% Convert datenum time to ticks (this is what my code spits out from Unity,
% plus it will be easier to fit a regression when the datetime is
% represented as a single value).
pre_time = nan(1,length(pre_time_orig));
post_time = nan(1,length(post_time_orig));
for i = 1:length(pre_time_orig)
    t = datenum(pre_time_orig{i});
    t = datenum2ticks(t);
    pre_time(i) = t;
end

for i = 1:length(post_time_orig)
    t = datenum(post_time_orig{i});
    t = datenum2ticks(t);
    post_time(i) = t;
end

% Plot pre-test pulses and post-test pulses
figure;
plot(pre_time,prePulseBins,'b*')
title('Pre-test pulses')

figure;
plot(post_time,postPulseBins,'r*')
title('Post-test pulses')

% % Calculate regressions for pre- and post-test timing
pre_P = polyfit(pre_time,prePulseBins',1);
post_P = polyfit(post_time,postPulseBins',1);

% Calculate regression using both 
all_time = cat(2,pre_time,post_time);
all_pulses = cat(1,prePulseBins,postPulseBins)';
all_P = polyfit(all_time,all_pulses,1);
all_y = polyval(all_P,all_time);
pre_y = polyval(pre_P,all_time);
post_y = polyval(post_P,all_time);

figure;
plot(all_time,all_pulses,'k*')
hold on;
plot(all_time,all_y,'r-')
title('All pulses')

% Save regression values for later
time_sync_regression = all_P;
save([subject_dir 'Mat Files/' save_stem  'time_sync.mat'],'time_sync_regression');