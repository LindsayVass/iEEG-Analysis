%UCDMC11 Preprocessing

%notes:
% 1) Data must be in ASCII format and must not have trigger information
% included when exported.  The trigger information adds a column to the
% data such that there are is one more column than channcels, throwing off
% the reshaping step.
% 2) Data must be in .edf format

clear all

%%% Set up paths
addpath(genpath('/Users/amberschedlbauer/Documents/MATLAB/EEG/iEEG/Subjects/UCDMC13/'))
%%%

%%% UDIs

subject_dir = '/Users/amberschedlbauer/Documents/MATLAB/EEG/iEEG/Subjects/UCDMC13/';
epoch_savefile = 'epoch_UCDMC13';

%cond1_labels = [111]; %cue spat correct
%cond2_labels = [121]; %cue temp correct

%%%

%%% Choose behavioral data file
cd(subject_dir);

cd Behavioral_data/Retrieval_data_uncorr/UCDMC13_112114/

% Behavioral file
d1 = dir;
cd_files = {d1.name};
[s,v] = listdlg('PromptString','Select Subject#_Data.mat file','SelectionMode','single','ListString',cd_files);
behav_filename = [subject_dir,'Behavioral_data/Retrieval_data_uncorr/UCDMC13_112114/',cd_files{s}]; % full path and file name

%%%

% %convert ascii to mat file, from T. Ellmore
% if exist([outfname,'.mat']) ~=2
%     mult = 10;  %factor to multiply and divide by, see convert_nkascii2mat for explanation from T. Ellmore
%     convert_nkascii2mat(fname,outfname,mult);
% else
%     load([outfname,'.mat'])
% end
% 
% %load converted mat file
% load(outfname)

% recode behavior because of wrong response keys used initially
cd(subject_dir);
load(behav_filename)
%[spat_acc,temporal_acc] = recode_houston_retrieval_behavior(behav_filename,'L','R');

%load patient LFP data in .edf format
cd Raw_data/
[hdr,data] = edfread('UCDMC13_112014.edf');

% %  From T. Ellmore convert_nkascii2mat documentation
% %  FYI: Data in the NK ASCII file are stored with floating point precision to the hundreth
% %  decimal place. On conversion, each data value is multiplied by 10 in order to store
% %  data as 16 bit integers (storage as 16 bit int rather than 32 bit float saves space).
% %  Therfore, when reading the resultant mat file each data value must be divided by the
% %  nkdata.multiplier, which is given by the command line parameter mult.
% %  You must do this to preserve micovolt units. Here's how:
% 
% %as recommended from above, divide by multiplier
% nkdata.eeg = nkdata.eeg./nkdata.multiplier;
% 
% %must reshape the data, as it is incorrectly shaped from the conversion
% %step above.
% 
% all_data = nkdata.eeg(:);
% 
% for n = 1:nkdata.nchannels
%     reshape_data(:,n) = nkdata.eeg(n,:);
% end

% Put in EEG struct
eeglab

EEG.setname = 'UCDMC13';
EEG.filename = 'UCDMC13';
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
hdr.label = cellfun(@(x) x(4:end-3),hdr.label(1:end-1),'UniformOutput',false);
hdr.label{hdr.ns} = 'marker';

% Work on channel locations, inserting dummy information so we can use ica
% and plot channel properties
ctr=1;

for e = 1:size(data,1)
    
    if isempty(strfind(hdr.label{e},'empty')); %only use real channels
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
eeglab redraw %updates changes in the EEG structure in the GUI

clear nkdata reshape_data

pop_eegplot(EEG,1,1,1);
keyboard
save('shortcut','EEG');

if exist('retrieval_onsets.mat','file') ~=2
    %request pulse channel #
    chz = {EEG.chanlocs.labels};
    [pulse_channel_number] = listdlg('PromptString','Select pulse channel:',...
        'SelectionMode','single','ListString',chz);
    
    %extract pulses. newer pulse extraction utilizing info from viewing data
    %in eeglab
    [realonset,interpolated] = Extract_sync_event_Houston_v3(EEG,pulse_channel_number,master_event_list,2);
    
    save('retrieval_onsets','realonset','interpolated','pulse_channel_number');
    
else
    load('retrieval_onsets.mat')
end

% Make vector of accuracy by duplicating and concatenating temporal and
% spatial_acc
full_spat_acc_vector = repmat(spat_acc,2,1);
full_spat_acc_vector = full_spat_acc_vector(:)';

full_temp_acc_vector = repmat(temporal_acc,2,1);
full_temp_acc_vector = full_temp_acc_vector(:)';

if strcmp(master_event_list.block{1},'Space')
    all_acc = [full_spat_acc_vector full_temp_acc_vector];
else
    all_acc = [full_temp_acc_vector full_spat_acc_vector];
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%% Conjunctive Coding
% Get proportion of triads in which subject got both spat and temp compared
% to only spat or only temp

for n = 1:3
    
    counter = 1;
    
    if n ==1
        store_data = master_event_list.refstore;
    elseif n ==2
        store_data = master_event_list.store1store;
    elseif n==3
        store_data = master_event_list.store2store;
    end
    
    for t = 1:2:140
        
        store = store_data(t);
        
        if strcmp(store,'candle_shop.png');
            recoded_stores_to_numbers(n,counter) = 1;
        elseif strcmp(store,'flower_patch_mod.png');
            recoded_stores_to_numbers(n,counter) = 2;
        elseif strcmp(store,'cookie_shop.png');
            recoded_stores_to_numbers(n,counter) = 3;
        elseif strcmp(store,'coffee_store.png');
            recoded_stores_to_numbers(n,counter) = 4;
        elseif strcmp(store,'family_place.png');
            recoded_stores_to_numbers(n,counter) = 5;
        end
        
        %block type vector
        if strcmp(master_event_list.block(t),'Time')
            recoded_stores_to_numbers(4,counter) = 1;
        elseif strcmp(master_event_list.block(t),'Space')
            recoded_stores_to_numbers(4,counter) = 2;
        end
        
        %accuracy
        recoded_stores_to_numbers(5,counter) = all_acc(t);
        
        counter = counter + 1;
        
    end
end

% Find each trial's triad buddy, independent of left-right status

for n = 1:size(recoded_stores_to_numbers,2)
    clear ref_idx store1_idx store2_idx tmp_intersect exact_final_intersect
    %find exact triad
    ref_idx = find(recoded_stores_to_numbers(1,:) == recoded_stores_to_numbers(1,n));
    store1_idx = find(recoded_stores_to_numbers(2,:) == recoded_stores_to_numbers(2,n));
    store2_idx = find(recoded_stores_to_numbers(3,:) == recoded_stores_to_numbers(3,n));
    
    tmp_intersect = intersect(ref_idx,store1_idx);
    exact_final_intersect = intersect(tmp_intersect,store2_idx);
    
    %find left-right switched triad
    clear ref_idx store1_idx store2_idx tmp_intersect switched_final_intersect
    
    ref_idx = find(recoded_stores_to_numbers(1,:) == recoded_stores_to_numbers(1,n));
    store1_idx = find(recoded_stores_to_numbers(3,:) == recoded_stores_to_numbers(2,n));
    store2_idx = find(recoded_stores_to_numbers(2,:) == recoded_stores_to_numbers(3,n));
    
    tmp_intersect = intersect(ref_idx,store1_idx);
    switched_final_intersect = intersect(tmp_intersect,store2_idx);

    all_matching_triads = [exact_final_intersect switched_final_intersect];
    
    %remove current trial, which will always match
    all_matching_triads = setdiff(all_matching_triads,n);
    
    %compare triads where one is space and the other is time
    opposite_triads = all_matching_triads(find(recoded_stores_to_numbers(4,all_matching_triads) ~= recoded_stores_to_numbers(4,n))); %true conjunctive

    if recoded_stores_to_numbers(5,n)==1 && sum(recoded_stores_to_numbers(5,opposite_triads)) >0 %%length(opposite_triads)   %at least one space and time
        trial_type_breakdown(n) = 1;
    elseif recoded_stores_to_numbers(5,n)==1 && sum(recoded_stores_to_numbers(5,opposite_triads)) ==0 && recoded_stores_to_numbers(4,n)==1 %time only
        trial_type_breakdown(n) = 2;
    elseif recoded_stores_to_numbers(5,n)==1 && sum(recoded_stores_to_numbers(5,opposite_triads)) ==0 && recoded_stores_to_numbers(4,n)==2 %space only
        trial_type_breakdown(n) = 3;
    end
    
end

fooD = repmat(trial_type_breakdown,2,1);
full_event_trial_type_breakdown = fooD(:)';


% Insert trigger events into data. Also, get indices for different conditions
event_counter = 1;
cue_spat_corr_idx = [];
cue_spat_incorr_idx = [];
cue_temp_corr_idx = [];
cue_temp_incorr_idx = [];
resp_spat_corr_idx = [];
resp_spat_incorr_idx = [];
resp_temp_corr_idx = [];
resp_temp_incorr_idx = [];

for n = 1:length(realonset)
    
    EEG.event(event_counter).latency = realonset(n);
    
    if strcmp(master_event_list.eventtype{event_counter},'Stores_onset') && strcmp(master_event_list.block{event_counter},'Space') && all_acc(event_counter) ==1 %cue spatial correct
        EEG.event(event_counter).type = '111';
        cue_spat_corr_idx = [cue_spat_corr_idx;event_counter];
    elseif strcmp(master_event_list.eventtype{event_counter},'Stores_onset') && strcmp(master_event_list.block{event_counter},'Space') && all_acc(event_counter) ==0 %cue spatial incorrect
        EEG.event(event_counter).type = '112';
        cue_spat_incorr_idx = [cue_spat_incorr_idx;event_counter];    
    elseif strcmp(master_event_list.eventtype{event_counter},'Stores_onset') && strcmp(master_event_list.block{event_counter},'Time') && all_acc(event_counter) ==1 %cue temporal correct
        EEG.event(event_counter).type = '121';
        cue_temp_corr_idx = [cue_temp_corr_idx;event_counter];
    elseif strcmp(master_event_list.eventtype{event_counter},'Stores_onset') && strcmp(master_event_list.block{event_counter},'Time') && all_acc(event_counter) ==0 %cue temporal incorrect
        EEG.event(event_counter).type = '122';
        cue_temp_incorr_idx = [cue_temp_incorr_idx;event_counter];
    elseif strcmp(master_event_list.eventtype{event_counter},'Response') && strcmp(master_event_list.block{event_counter},'Space') && all_acc(event_counter) ==1 %response spatial correct
        EEG.event(event_counter).type = '211';
        resp_spat_corr_idx = [resp_spat_corr_idx;event_counter];
    elseif strcmp(master_event_list.eventtype{event_counter},'Response') && strcmp(master_event_list.block{event_counter},'Space') && all_acc(event_counter) ==0 %response spatial incorrect
        EEG.event(event_counter).type = '212';
        resp_spat_incorr_idx = [resp_spat_incorr_idx;event_counter];
    elseif strcmp(master_event_list.eventtype{event_counter},'Response') && strcmp(master_event_list.block{event_counter},'Time') && all_acc(event_counter) ==1 %response temporal correct
        EEG.event(event_counter).type = '221';
        resp_temp_corr_idx = [resp_temp_corr_idx;event_counter];
    elseif strcmp(master_event_list.eventtype{event_counter},'Response') && strcmp(master_event_list.block{event_counter},'Time') && all_acc(event_counter) ==0 %response temporal incorrect
        EEG.event(event_counter).type = '222';
        resp_temp_incorr_idx = [resp_temp_incorr_idx;event_counter];
    end
   
    event_counter = event_counter + 1;
end

num_conjunctive = length(find(trial_type_breakdown==1));
num_time = length(find(trial_type_breakdown==2));
num_space = length(find(trial_type_breakdown==3));

% No localization stuff so not doing it yet

% Re-reference the data
for iElec = 1:size(EEG.data,1)-1 %no marker
    EEG.data(iElec,:) = EEG.data(iElec,:)-EEG.ref_tseries;
end
eeglab redraw

% Save unepoched data from entire file with all triggers but otherwise fully processed
unepoch_savefile = ['un',epoch_savefile];

if exist([unepoch_savefile,'.set']) ~=2
    EEG = pop_saveset(EEG, 'filename', unepoch_savefile)%,'filepath',str2);
    eeglab redraw
end

% Make  cue epochs
[EEG] = pop_epoch(EEG,{},[-1 2.2]);

% Save data
EEG = pop_saveset(EEG, 'filename', epoch_savefile);
eeglab redraw