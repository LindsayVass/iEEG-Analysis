function PutIntoEEG(subject_dir, subject_id, edf_files,thisEDF, save_stems)
% Set up EEG structure from EDF data. Will save a .mat file containing the
% EEG data with bad channels removed.
%
% PutIntoEEG(subject_dir, subject_id, edf_files, thisEDF, save_stems)
%
% INPUTS:
%   subject_dir: main directory for subject
%   subject_id: ID for subject (e.g., UCDMC15)
%   edf_files: cell array containing paths to the EDF files
%   thisEDF: which EDF we're pre-processing
%   save_stems: cell array of stems for saving
%

% Load patient LFP data in .edf format
[hdr,data] = edfread(edf_files{thisEDF});

%% Put in EEG struct
eeglab

EEG = eeg_emptyset();

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

% Work on channel locations, inserting dummy information
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

%% Find marker channel and remove bad channels


% Prompt user to select marker channel number
input('\n\n\nIdentify the marker channel and press ENTER when ready to continue.')
prompt = 'Select marker channel';
[markerChanNumber, ~] = listdlg('ListString',{EEG.chanlocs.labels},'PromptString',prompt);


markerChanName = hdr.label(markerChanNumber);
EEG.chanlocs(markerChanNumber).labels = 'marker';



eeglab redraw % updates changes in the EEG structure in the GUI
pop_eegplot(EEG,1,1,1);
% Prompt user for which channels to remove
input('\n\n\nIdentify the bad channels to remove and press ENTER when ready to continue.')
prompt = 'Select bad channels to remove';
[badChans, ~] = listdlg('ListString',{EEG.chanlocs.labels},'PromptString',prompt);


badChanNames = {EEG.chanlocs([badChans]).labels};
EEG.data(badChans,:) = [];
EEG.nbchan = EEG.nbchan - length(badChans);
EEG.chanlocs([badChans]) = [];



eeglab redraw % updates changes in the EEG structure in the GUI
pop_eegplot(EEG,1,1,1);
% Prompt user to identify channels showing spiking activity
input(['\n\n\nIdentify the channels showing the most prominent spiking activity.' ...
    '\nWe will hang on to them for now, but they will not be included in the' ...
    '\naverage when we re-reference the data. Press ENTER when ready to continue.'])
prompt = 'Select channels showing prominent spiking activity.';
[spikingChans, ~] = listdlg('ListString',{EEG.chanlocs.labels},'PromptString',prompt);


spikingChanNames = {EEG.chanlocs([spikingChans]).labels};

% Save current state of EEG
save([subject_dir 'PreProcessing Intermediates/' save_stems{thisEDF} '_badChansRemoved.mat'],'EEG','markerChanName','badChanNames','markerChanNumber','badChans','spikingChans','spikingChanNames');

end