function [systemTime, xPos, zPos] = parseBehavioralTxt(inputPath, fixFile)
% [systemTime, xPos, zPos] = parseBehavioralTxt(inputPath)
% [systemTime, xPos, zPos] = parseBehavioralTxt(inputPath, fixFile)
% 
% Purpose: Read in the txt file output by Unity during the experiment and
%   extract the columns corresponding to the system time in ticks, as well
%   as the avatar's X and Z coordinates
%
% INPUT:
%   inputPath: path to the txt file
%
% OPTIONAL INPUT:
%   fixFile: set to 1 if this file was manually fixed (this is true for
%       UCDMC15)
%
% OUTPUT:
%   systemTime: vector of values representing the system time in Unix ticks
%   xPos: vector of values representing the avatar's X coordinate
%   zPos: vector of values representing the avatar's Z coordinate
%
% Author: Lindsay Vass
% Date: 5 August 2015

%% Load in the unity output
fid  = fopen(inputPath);
if ~exist('fixFile', 'var')
    data =  textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
else
    data = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n');
end
fclose(fid);

%% Extract variables of interest
systemTime   = data{3};
xPos         = data{8};
zPos         = data{9};

end