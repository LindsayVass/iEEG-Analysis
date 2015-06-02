% Build data structure for analysis
%
% Lindsay Vass
% 29 May 2015


field1 = 'subjectID';
field2 = 'teleporter';
field3 = 'numEDFs';
field4 = 'chanList'; 
sessionInfo = struct(field1, {}, field2, {}, field3, {}, field4, {});


% UCDMC13
sessionInfo(1).subjectID  = 'UCDMC13';
sessionInfo(1).teleporter = {'TeleporterA'};
sessionInfo(1).numEDFs    = [1];
sessionInfo(1).chanList   = {'LAD1' 'LHD1'};

% UCDMC14
sessionInfo(2).subjectID  = 'UCDMC14';
sessionInfo(2).teleporter = {'TeleporterA', 'TeleporterB'};
sessionInfo(2).numEDFs    = [1 1];
sessionInfo(2).chanList   = {'LAD1' 'LHD1' 'LHD2' 'RAD1' 'RHD1' 'RHD2'};

% UCDMC15
sessionInfo(3).subjectID  = 'UCDMC15';
sessionInfo(3).teleporter = {'TeleporterA', 'TeleporterB'};
sessionInfo(3).numEDFs    = [2 2];
sessionInfo(3).chanList   = {'RAD3' 'RAD4' 'RAD5' 'RAD6' 'RHD1' 'RHD2' 'RHD3' 'RHD4' 'LAD3' 'LAD4' 'LAD5' 'LHD1' 'LHD2' 'LHD3'};

save('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo.mat', 'sessionInfo');
