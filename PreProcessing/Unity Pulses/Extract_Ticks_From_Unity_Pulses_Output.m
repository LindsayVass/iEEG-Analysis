% This script will read in the pulses text file output from Unity and
% extract the time in ticks for each pulse onset.
%
% Lindsay Vass 24 April 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set file path
fileName = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Behavioral Data/TeleporterA/s3_PULSES_.txt';
saveFile = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Raw Data/pulses_teleporterA.mat';

%% Parse the file

fid = fopen(fileName);
pulses = [];

% Read in the pulses file line by line
while 1
    % Get the next line
    tline = fgetl(fid);
    
    % Check if this is a pulse start
    ind = strfind(tline,'Pulse Start');
    
    if ~isempty(ind) % If this is a pulse
        pulseText = tline(ind+12:end);
        pulses(end+1,1) = str2num(pulseText);
    end
    
    
    if ~ischar(tline)
        break
    end
    
    disp(tline)
end

fclose(fid);

% Ask if last pulse is partial
prompt = 'Does last pulse start but not end? Press Y to remove it: ';
answer = input(prompt,'s');

if(strcmpi('Y',answer) == 1)
    
    disp('Last pulse removed.');
    
    % Last pulse is partial, so remove it
    pulses(end) = [];
end

figure;
plot(pulses);
title('Pulse timing');

save(saveFile,'pulses');
