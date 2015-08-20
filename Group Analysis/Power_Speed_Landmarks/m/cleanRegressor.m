function [regressor1, regressor2] = cleanRegressor(inputRegressor, EDF1inds, EDF2inds, goodEpochs1, goodEpochs2, removedOnsets1, removedOnsets2)
% [regressor1, regressor2] = cleanRegressor(inputRegressor, EDF1inds, EDF2inds, goodEpochs1, goodEpochs2, removedOnsets1, removedOnsets2)
%
% Purpose: Take a vector of values and split it into two according to which
%   EDF the values come from. Then, keep only the valid (non-artifact)
%   values from each list.
%
% INPUT:
%   inputRegressor: vector of values to clean; could be either numeric
%       (vector) or string (cell array)
%   EDF1inds: vector of integers indicating the indices that belong to EDF1
%   EDF2inds: vector of integers indicating the indices that belong to EDF2
%   goodEpochs1: vector of integers indicating which EDF1inds are valid
%   goodEpochs2: vector of integers indicating which EDF2inds are valid;
%       for example, if EDF2inds is [25 26 27] and goodEpochs2 is [2], it 
%       indicates that "26" is valid whereas "25" and "27" are not
%   removedOnsets1: indices of EDF1 epochs that were removed prior to
%       calculating goodEpochs because they overlap with teleportation
%   removedOnsets2: indices of EDF2 epochs that were removed prior to
%       calculating goodEpochs because they overlap with teleportation 
%
% OUTPUT:
%   regressor1: vector or cell array of values for valid trials in EDF1
%   regressor2: vector or cell array of values for valid trials in EDF2
%
% Author: Lindsay Vass
% Date: 5 August 2015

regressor1 = inputRegressor(EDF1inds);
regressor2 = inputRegressor(EDF2inds);

regressor1(removedOnsets1) = [];
regressor2(removedOnsets2) = [];

regressor1 = regressor1(goodEpochs1);
regressor2 = regressor2(goodEpochs2);

end