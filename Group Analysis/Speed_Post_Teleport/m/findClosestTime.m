function [ind, timeDiff] = findClosestTime(testTime, timeList)
% findClosestTime: given a vector of times, find the index of the time
% in the list closest to the time of interest
% >> [ind, timeDiff] = findClosestTime(testTime, timeList)
%
% INPUTS:
%   testTime: value of interest
%   timeList: vector of values
%
% OUTPUTS:
%   ind: index of the closest value in timeList
%   timeDiff: difference between testTime and timeList(ind)

absTimeDiff = abs(timeList - testTime);
ind         = find(absTimeDiff == min(absTimeDiff));
timeDiff    = testTime - timeList(ind);

end
