% set file path
fileName = '/Users/Lindsay/Documents/MATLAB/iEEG/Subjects/UCDMC14/Behavioral Data/s2_patientTeleporterData 2.txt';

%% Parse the file

[~,timeElapsed,systemTime,~,~,~,~,~,~,~] = textread(fileName,'%n%n%n%s%s%s%s%n%n%n','delimiter',',','headerlines',1);


% Calculate tick interval
systemTime_shift = systemTime;
systemTime_shift(end) = [];
systemTime_shift = cat(1,nan,systemTime_shift);
tick_interval = systemTime - systemTime_shift;

% Look at the data to see if there are any weird discontinuities
hist(tick_interval);
