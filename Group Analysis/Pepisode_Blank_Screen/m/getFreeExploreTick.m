function freeExploreTick = getFreeExploreTick(freeExplorePath, isUCDMC15)
% Extract the last time point (in ticks) from the Free Explore text file.
% >> freeExploreTick = getFreeExploreTick(freeExplorePath, isUCDMC15)
%
% Inputs:
%   freeExplorePath: path to the text file
%   isUCDMC15: boolean for whether this subject is UCDMC15 (who has a
%       different file format)
%
% Output:
%   freeExploreTick: final time point of the text file in Unity ticks



% Load in the unity output
fid  = fopen(freeExplorePath);
if isUCDMC15 == 0
    data =  textscan(fid,'%f%f%f%f%f','delimiter',',','Headerlines',1);
else
    data = textscan(fid,'%f%f%f%f%f','delimiter',',','Headerlines',1, 'EndOfLine','\r\n');
end

fclose(fid);

% Extract variable of interest
systemTime = data{2};

% Extract final time point
freeExploreTick = systemTime(end);
