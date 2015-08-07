function [newRegressor] = factorRegressor(stringRegressor, mappings)
% [newRegressor] = factorRegressor(stringRegressor, mappings)
%
% Purpose: Take a categorical regressor with multiple levels and convert
%   them to integers so that they can be entered into a regression model.
%
% INPUT:
%   stringRegressor: cell array of strings
%   mappings: structure with two fields, "name" which indicates the strings
%       found in stringRegressor, and "value" which indicates the integers
%       the strings should be mapped to
%
%   Example mappings:
%       mappings = struct('name', {'RICH', 'POOR', 'CENTRAL'}, ...
%                         'value', {2, 1, 2});
%
% OUTPUT:
%   newRegressor: vector of integers representing the values in
%       stringRegressor
%
% Author: Lindsay Vass
% Date: 6 August 2015

newRegressor = nan(size(stringRegressor));
for thisFactor = 1:size(mappings, 2)
    thisName  = mappings(thisFactor).name;
    thisValue = mappings(thisFactor).value;
    factorInd = strcmpi(thisName, stringRegressor);
    newRegressor(factorInd) = thisValue;
end

end