function [realonsets,interpolated] = Extract_sync_event_Houston_v3(EEG,pulse_channel_number,master_event_list,viewer_type)

% Identical to v2 but allows for inspecting the data in EEGlab versus
% nkviewer, since the more recent data is already exported by Chris

%fixed bug with use_idx, where selecting a index returned the index of the
%index, rather than the selected number.  

%EEG:  EEG data.  Must include the field EEG.start_time
%pulse_channel_number
% viewer_type is a switch between nkviewer (1) and eeglab (2)

tries_counter = 1;
retry = 1;
interpolated(1) = 0;  %first pulse cannot be interpolated

% %convert timeshift to samples
% use_timeshift = round(timeshift*EEG.srate); %samples

%make EEG_Start a time/date vector which can be used to determine elapsed
first_logged_pulse = master_event_list.time{1};
EEG_Start = first_logged_pulse;

%modified to handle case where hour information is a single number.  
%AJW  4-24-13
% if isnumeric(EEG.start_time(2)) 
% EEG_Start(4) = str2num(EEG.start_time(1:2));
% EEG_Start(5) = str2num(EEG.start_time(4:5));
% EEG_Start(6) = str2num(EEG.start_time(7:8));
% else
% EEG_Start(4) = str2num(EEG.start_time(1));
% EEG_Start(5) = str2num(EEG.start_time(3:4));
% EEG_Start(6) = str2num(EEG.start_time(6:7));
% end

%%% Amber and Millie
[hour,minute] = strtok(EEG.start_time,':');
minute = minute(2:end);
[minute,sec] = strtok(minute,'.'); %IMPORTANT: SOMETIMES ITS A ':' OR A '.'. LOOK!!!
sec = sec(2:end);

EEG_Start(4) = str2num(hour);
EEG_Start(5) = str2num(minute);
EEG_Start(6) = str2num(sec);
% if numel(sec) == 1
%     EEG_Start(6) = ((str2num(sec)*60))/(100);
% else
%     EEG_Start(6) = ((str2num(sec)*60))/(10^(numel(sec)));
% end
%%%

% Request UI for first retrieval pulse based on visual inspection of data
if viewer_type ==1 %nkdata with hour, min, sec
    prompt = {'NK Pulse hour','NK Pulse Min','NK Pulse Sec'};
    dlg_title = 'Inputs from inspecting data NK Viewer';
    num_lines = 1;
    def = {'0','0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    hr = str2num(cell2mat(answer(1)));
    after_buffer = str2num(cell2mat(answer(2)));
    
    first_viewed_pulse = first_logged_pulse;
    first_viewed_pulse(4) = str2num(cell2mat(answer(1)));
    first_viewed_pulse(5) = str2num(cell2mat(answer(2)));
    first_viewed_pulse(6) = str2num(cell2mat(answer(3)));
else
    
    prompt = {'Add hours','Add Mins','Add Sec'};
    dlg_title = 'Inputs from inspecting data in EEGLab';
    num_lines = 1;
    def = {'0','0','0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    add_seconds = str2num(cell2mat(answer(1)));
    
    first_viewed_pulse = EEG_Start;
    first_viewed_pulse(4) = EEG_Start(4)+str2num(cell2mat(answer(1)));
    first_viewed_pulse(5) = EEG_Start(5)+str2num(cell2mat(answer(2)));
    first_viewed_pulse(6) = EEG_Start(6)+str2num(cell2mat(answer(3)));
    
end

init_elapse = etime(first_viewed_pulse,EEG_Start); %Computes the time in seconds between two date vectors.

%get initial time windows
%init_elapse = etime(master_event_list.time{1},EEG_Start)+timeshift; %incorporate timeshift in
init_indices = -EEG.srate+round(init_elapse*EEG.srate):EEG.srate+round(init_elapse*EEG.srate);

figure;plot(EEG.data(pulse_channel_number,init_indices))

%disp(['First trigger is around ',num2str(init_elapse),' seconds'])

%request pulse polarity
[pulse_polarity,v] = listdlg('PromptString','Pulse polarity?','SelectionMode','single','ListString',{'Positive','Negative'});

%Get initial parameters
prompt = {'Before trim amount:','After trim amount:'};
dlg_title = 'Input';
num_lines = 1;
def = {num2str(-EEG.srate),num2str(EEG.srate)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
before_buffer = str2num(cell2mat(answer(1)));
after_buffer = str2num(cell2mat(answer(2)));

while retry == 1
    
    if tries_counter > 1
        
        %request UI for first retrieval pulse based on visual inspection of
        prompt = {'NK Pulse hour','NK Pulse Min','NK Pulse Sec'};
        dlg_title = 'Inputs from inspecting data NK Viewer';
        num_lines = 1;
        def = {'0','0','0'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        hr = str2num(cell2mat(answer(1)));
        after_buffer = str2num(cell2mat(answer(2)));
        
        first_viewed_pulse = [];
        first_viewed_pulse = first_logged_pulse;
        first_viewed_pulse(4) = str2num(cell2mat(answer(1)));
        first_viewed_pulse(5) = str2num(cell2mat(answer(2)));
        first_viewed_pulse(6) = str2num(cell2mat(answer(3)));
        
        prompt = {'Before trim amount:','After trim amount:'};
        dlg_title = 'Input';
        num_lines = 1;
        def = {num2str(before_buffer),num2str(after_buffer)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        before_buffer = str2num(cell2mat(answer(1)));
        after_buffer = str2num(cell2mat(answer(2)));
    end
    
    index_map = [];
    indices_to_check = [];
    tmp_data = [];
    for pulz = 1:length(master_event_list.time)
        
        %get amount of time elapsed into recording for each pulz
        if pulz ==1
            elapsed(pulz) = init_elapse;
        else
            first_pulse_idx = realonsets(1);
            elapsed(pulz) = init_elapse+etime(master_event_list.time{pulz},master_event_list.time{1});
        end
        
        %determine indice window
        
        indices_to_check(pulz,:) = before_buffer+round(elapsed(pulz)*EEG.srate):after_buffer+round(elapsed(pulz)*EEG.srate);
        %get pulse data
        try
        tmp_data(pulz,:) = EEG.data(pulse_channel_number,indices_to_check(pulz,:));
        catch
            keyboard
        end
        
        if max(abs(tmp_data(pulz,:))) > 50  %look for a pulse, positive or negative
            interpolated(pulz) = 0;  %didnt interpolate this trial
            if pulse_polarity ==1 %positive
                [idx ] = find(tmp_data(pulz,:)==max(tmp_data(pulz,:)));
            elseif pulse_polarity ==2 %negative
                [idx ] = find(tmp_data(pulz,:) == min(tmp_data(pulz,:)));
            else
                error('Pulses cant be zero you buffoon');
            end
        else  %if it can't find something , try to interpolate.  this step is done down below when realonsets is defined
            
            interpolated(pulz) = 1;
        end
        
        if length(idx) > 1 % if it finds more than one, do some stuff
            %plot image data so far
            figure;
            subplot(3,1,1);
            imagesc(tmp_data);
            title('Real data')
            subplot(3,1,2);
            image(index_map)
            title('Extracted indices')
            
            %plot EEG data around time point
            subplot(3,1,3);
            plot(EEG.data(pulse_channel_number,indices_to_check(pulz,:)));
            %prompt for choice
            
            %
            strngs = {num2str(idx'),'None'};
            [use_idx] = listdlg('PromptString','Select pulse onset:',...
                'SelectionMode','single','ListString',strngs);
            
            % Amber 5/19/14: Moved from line 200 
            if use_idx == length(idx)+1; %implies NONE was selected
                break
            end
            %
            
            use_idx = str2num(strngs{1}(use_idx,:));  %added 4/24/2013.  use_idx above is the index of the values selected, rather than the actual number.
            %so if you have the options 30,31,32,33 and you pick "30",
            %use_idx = 1, not 30, which is a problem
            
        else %otherwise use what it found
            use_idx = idx;
        end
        
        %specify indices used and realonsets
        if interpolated(pulz) ==0  %don't interpolate, use what it found
            realonsets(pulz) = indices_to_check(pulz,use_idx);
            index_map(pulz,:) = zeros(length(tmp_data(pulz,:)),1);
            index_map(pulz,use_idx) = 100;
        
        else  %interpolate
            elapsed_time = etime(master_event_list.time{pulz},master_event_list.time{1});
            realonsets(pulz) = first_pulse_idx+round(elapsed_time*EEG.srate);
            index_map(pulz,:) = ones(length(tmp_data(pulz,:)),1).*100;  %all 100s to show interpolated
        
        end
        close all
    end
    
    %plot final data and indices extracted
    figure;
    subplot(2,1,1);
    imagesc(tmp_data);
    title('Real data')
    subplot(2,1,2);
    image(index_map)
    title('Extracted indices')
    
    %prompt for rebate
    prompt = {'Retry? (1=yes,0=no)'};
    dlg_title = 'Input';
    num_lines = 1;
    def = {'1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    retry = str2num(cell2mat(answer(1)));
    tries_counter =tries_counter + 1;
end
close all
end
