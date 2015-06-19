% inputPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Behavioral Data/TeleporterB/s3_FindStore_TeleporterB.txt';
% outputPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Behavioral Data/TeleporterB/s3_FindStore_TeleporterB_FIXED.txt';

inputPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Behavioral Data/TeleporterA/s3_FindStore_TeleporterA.txt';
outputPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC15/Behavioral Data/TeleporterA/s3_FindStore_TeleporterA_FIXED.txt';

fid  = fopen(inputPath);

txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n');

fclose(fid);

trialNumber = txtdata{1};
timeElapsed = txtdata{2};
systemTime = txtdata{3};
target = txtdata{4};
armType = txtdata{5};
spaceType = txtdata{6};
timeType = txtdata{7};
xPos = txtdata{8};
zPos = txtdata{9};
yRot = txtdata{10};

for thisRow = 1:length(armType)
   
    switch target{thisRow}
        case 'Coffee Shop' 
            armType{thisRow} = 'Rich';
            spaceType{thisRow} = 'Far';
        case 'Ice Cream'
            armType{thisRow} = 'Rich';
            spaceType{thisRow} = 'Far';
        case 'Arm1 Teleporter'
            armType{thisRow} = 'Rich';
            spaceType{thisRow} = 'Far';
        case 'Florist'
            armType{thisRow} = 'Poor';
            spaceType{thisRow} = 'Near';
        case 'Music Store'
            armType{thisRow} = 'Poor';
            spaceType{thisRow} = 'Near';
        case 'Arm2 Teleporter'
            armType{thisRow} = 'Poor';
            spaceType{thisRow} = 'Near';
        case 'Grocery Store' 
            armType{thisRow} = 'Rich';
            spaceType{thisRow} = 'Near';
        case 'Clothing Store'
            armType{thisRow} = 'Rich';
            spaceType{thisRow} = 'Near';
        case 'Arm3 Teleporter'
            armType{thisRow} = 'Rich';
            spaceType{thisRow} = 'Near';
        case 'Pet Store'
            armType{thisRow} = 'Poor';
            armType{thisRow} = 'Far';
        case 'Dentist'
            armType{thisRow} = 'Poor';
            armType{thisRow} = 'Far';
        case 'Arm4 Teleporter'
            armType{thisRow} = 'Poor';
            armType{thisRow} = 'Far';
            
    end
    
end


fid = fopen(outputPath, 'w');
fprintf(fid, '%s\n', 'TrialNumber,TimeElapsed,SystemTime,Target,ArmType,TeleporterSpaceType,TeleporterTimeType,XPosition,ZPosition,YRotation');

for thisRow = 1:length(trialNumber)
    
    fprintf(fid, '%d,%f,%18.0f,%s,%s,%s,%s,%3.4f,%3.4f,%1.10f\n',trialNumber(thisRow),timeElapsed(thisRow),systemTime(thisRow),target{thisRow},armType{thisRow},spaceType{thisRow},timeType{thisRow},xPos(thisRow),zPos(thisRow),yRot(thisRow));
    
end

fclose(fid)