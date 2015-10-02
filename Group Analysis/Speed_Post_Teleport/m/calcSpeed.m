function meanSpeed = calcSpeed(xVal, zVal, txtIntervalMs)

distMat = dist([cat(1, xVal', zVal')]);

% restrict to distances between consecutive points
p1 = [1:1:length(xVal)-1];
p2 = [2:1:length(xVal)];
distances = distMat(sub2ind(size(distMat), p1, p2));
meanSpeed = mean(distances / txtIntervalMs);

end