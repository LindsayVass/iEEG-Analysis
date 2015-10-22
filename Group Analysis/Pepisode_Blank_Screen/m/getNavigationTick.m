function navigationTick = getNavigationTick(navigationPath, isUCDMC15)
% Extract the last time point (in ticks) from the Navigation text file.
% >> navigationTick = getNavigationTick(navigationPath, isUCDMC15)
%
% Inputs:
%   navigationPath: path to the text file
%   isUCDMC15: boolean for whether this subject is UCDMC15 (who has a
%       different file format)
%
% Output:
%   navigationTick: final time point of the text file in Unity ticks



% Load in the unity output
fid  = fopen(navigationPath);
if isUCDMC15 == 0
    data =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
else
    data = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1, 'EndOfLine','\r\n');
end

fclose(fid);

% Extract variable of interest
systemTime = data{3};

% Extract final time point
navigationTick = systemTime(end);
