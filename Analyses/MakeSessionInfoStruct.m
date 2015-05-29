% Build data structure for analysis
%
% Lindsay Vass
% 29 May 2015


field1 = 'subjectID';     value1 = {'UCDMC13','UCDMC14','UCDMC15'};
field2 = 'teleporter';    value2 = {{'TeleporterA'}, ...
                                    {'TeleporterA', 'TeleporterB'}, ...
                                    {'TeleporterA', 'TeleporterB'}};
field3 = 'chanList';      value3 = {{'LAD1' 'LHD1'}, ...
                                    {'LAD1' 'LHD1' 'LHD2' 'RAD1' 'RHD1' 'RHD2'}, ...
                                    {'RAD3' 'RAD4' 'RAD5' 'RAD6' 'RHD1' 'RHD2' 'RHD3' 'RHD4' 'LAD3' 'LAD4' 'LAD5' 'LHD1' 'LHD2' 'LHD3'}};

sessionInfo = struct(field1, value1, field2, value2, field3, value3);

save('/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Subject Info/SessionInfo.mat', 'sessionInfo');
