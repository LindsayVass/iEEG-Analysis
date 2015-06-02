function chanNames = findChansOnThisElectrode(chanList, depthName)
% chanNames = findChansOnThisElectrode(chanList, depthName)
%
% Given a list of channel names (e.g., {'LAD2','LAD3','RAD2','RAD3'}),
% return the list of channel names that are on the current electrode of
% interest (e.g., 'LAD').
%
% INPUTS:
%   chanList: cell array of channel names; the first three letters of the
%       channel name = depthName
%   depthName: string corresponding to name of the depth electrode
%
% OUTPUT:
%   chanNames: cell array of channels from chanList that are on the depth
%       electrode named depthName
%

chanDepth = char(chanList);
chanDepth = chanDepth(:, 1:3); % keep the first 3 characters, removing the electrode #
chanDepth = cellstr(chanDepth);

chanInd = strcmpi(depthName, chanDepth);
chanNames = chanList(chanInd);

end