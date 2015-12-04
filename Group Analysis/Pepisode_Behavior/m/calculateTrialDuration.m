function calculateTrialDuration(inputPath, fixFile, subjectID, teleporter, outputPath)

% load the unity text file
fid  = fopen(inputPath);

if fixFile == 1
    txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1,'EndOfLine','\r\n');
else
    txtdata = textscan(fid,'%d%f%f%s%s%s%s%f%f%f','delimiter',',','Headerlines',1);
end

fclose(fid);

trialNumber = txtdata{1};
timeElapsed = txtdata{2};

% get indices of new trials
newTrialInd = [1; find(diff(trialNumber)) + 1; length(trialNumber)];

% get times of new trials
trialStartTimes = timeElapsed(newTrialInd);

% get durations
trialDurations = diff(trialStartTimes);
trialDurCell = mat2cell(trialDurations, ones(length(trialDurations), 1), 1);

output = cell(length(trialDurations) + 1, 3);
output(1,:) = {'SubjectID', 'Session', 'Duration'};
output(2:end, 1) = {subjectID};
output(2:end, 2) = {teleporter};
output(2:end, 3) = trialDurCell;

% save the output
dlmcell(outputPath, output, 'delimiter', ',');