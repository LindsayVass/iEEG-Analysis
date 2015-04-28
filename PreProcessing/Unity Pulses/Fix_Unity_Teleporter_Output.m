% Unity script did not correctly output the "ArmType" or
% "TeleporterSpaceType" variables, so we're going to fix the files here.
% Since these two variables are dependent on the maze itself, we can just
% grab the correct output from a previous subject's data:
%
% Lindsay Vass 21 April 2015

%% Set up script
clc;clear;close all;

% Path to the file
subject_dir = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/';
dataPath = [subject_dir 'Behavioral Data/TeleporterA/s3_FindStore_TeleporterA.txt'];

% file to save
saveFile = [subject_dir 'Behavioral Data/TeleporterA/s3_FindStore_TeleporterA_FIXED.txt'];

% List of targets along with their associated arm/space types
targetList = {'Coffee Shop' 'Florist' 'Grocery Store' 'Pet Store' 'Arm1 Teleporter' 'Arm2 Teleporter' 'Arm3 Teleporter' 'Arm4 Teleporter'};
armTypeList = {'Rich' 'Poor' 'Rich' 'Poor' 'Rich' 'Poor' 'Rich' 'Poor'};
spaceTypeList = {'Far' 'Near' 'Near' 'Far' 'Far' 'Near' 'Near' 'Far'};

 %% Load data
 fid  = fopen(dataPath);
 data = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
 fclose(fid);
 
 trialNumber = data{1};
 timeElapsed = data{2};
 systemTime  = data{3};
 target      = data{4};
 armType     = data{5};
 spaceType   = data{6};
 timeType    = data{7};
 xPos        = data{8};
 zPos        = data{9};
 yRot        = data{10};
 
 % open file to save
 fid = fopen(saveFile,'w');
 fprintf(fid,'%s\n','TrialNumber,TimeElapsed,SystemTime,Target,ArmType,TeleporterSpaceType,TeleporterTimeType,XPosition,ZPosition,YRotation');

 fixedData = cell(length(trialNumber),10);
for thisLine = 1:length(armType)
    
    % Get the current target
    thisTarget = target{thisLine};
    
    % Find the target in the list
    ind = strcmpi(thisTarget,targetList);
    
    % Apply the arm and space types based on targetList ind
    armType(thisLine) = {armTypeList(ind)};
    spaceType(thisLine) = {spaceTypeList(ind)};
    
    fixedData(thisLine,:) = {num2str(trialNumber(thisLine)),num2str(timeElapsed(thisLine)),sprintf('%18.0f',systemTime(thisLine)),target{thisLine},armType{thisLine},spaceType{thisLine},timeType{thisLine},num2str(xPos(thisLine)),num2str(yRot(thisLine)),num2str(zPos(thisLine))};
    
end

dlmcell(saveFile,fixedData,',','-a');

fclose(fid);
