function [rayleighP, rayleighZ] = multRayleighTest(phaseData)
% multRayleighTest: at each time point, run the Rayleigh Test to determine
% if phase distribution is non-uniform
% >> [rayleighP, rayleighZ] = multRayleighTest(phaseData)
%
% Input:
%   phaseData: trials x timepoints matrix of phase data
%
% Outputs:
%   rayleighP: vector of P values for each timepoint
%   rayleighZ: vector of Z values for each timepoint
%
% Lindsay Vass
% 15 September 2015

rayleighP = nan(1, size(phaseData, 2));
rayleighZ = rayleighP;
for i = 1:size(phaseData, 2)
    [p, z] = circ_rtest(phaseData(:, i));
    rayleighP(i) = p;
    rayleighZ(i) = z;
end