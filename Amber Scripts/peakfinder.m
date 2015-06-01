% Instead of Drew's code, this is another way to find the peaks in the
% marker channel. It is a little buggy because sometimes it finds local
% maxima that are right next to each other. This is really for a rough
% estimate.

%%%%% FOR SUBJECT 439!!! %%%%%

EEG2 = EEG;

EEG2.data = double(EEG2.data);

EEG2 = pop_resample(EEG,250);

[pks,idx_data] = findpeaks(EEG2.data(121,:),'minpeakheight',200); %The min peak height changes for each data file.

% Finds repeat peaks
idx = idx_data/250;
repeats = {};
round_idx = floor(idx);
uni_idx = unique(round_idx);
n = histc(round_idx,uni_idx);
multiple = find(n > 1);
for i_rep = 1:length(multiple)
    repeats{i_rep} = find(round_idx == uni_idx(multiple(i_rep)));
end

to_delete = cellfun(@(x) x(2:end),repeats,'UniformOutput',false);
bad_idx = cell2mat(to_delete);
idx(bad_idx) = [];
idx = [idx(1:3) idx(3)+1.4337 idx(4)-1.4337 idx(4:5) idx(5)+1.5087 idx(5)+1.5087+5.5425 idx(6) 160.4057 idx(7:25) 283.7190 idx(26:28) 296.7822 298.0591 idx(29:33) 323.5330 idx(34:end)];
%peaks = EEG2.data(idx_data);
time_between = diff(idx);
%time_between2 = time_between(134:end); %Index where retrieval pulses should supposedly start
time_bet = [time_between(2:2:46) time_between(47:2:end)];
RT = cell2mat(master_event_list.RT);
figure
plot(time_bet,'r')
hold
plot(RT,'b')

