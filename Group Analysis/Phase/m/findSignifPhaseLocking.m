% This script will examine the results of "runPhaseLockingAnalysis.m" to
% determine whether there is evidence of significant sustained phase
% locking. This is defined as a significant P value on the Rayleigh test
% for all samples corresponding to 2 cycles of the test frequency, within
% the 500 ms sampled interval.
%
% Lindsay Vass
% 16 September 2015

%phaseResultsPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Phase/mat/phaseResults.mat';
phaseResultsPath = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Phase/mat/phaseResults_rayleighP01.mat';

%savePath         = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Phase/csv/signifPhaseResults.csv';
savePath         = '/Users/Lindsay/Documents/MATLAB/iEEG/Group Analysis/Phase/csv/signifPhaseResults_rayleighP01.csv';

load(phaseResultsPath);
samplingRate   = params.samplingRate;
rayleighCycles = params.rayleighCycles;

outputCell = cell(1, 7);
outputCell(1, :) = {'Subject', 'Teleporter', 'Electrode', 'Analysis', 'Frequency', 'Onset', 'Offset'};
for thisSubject = 1:length(phaseResults)
    for thisTele = 1:length(phaseResults(thisSubject).teleporter)
        for thisDepth = 1:length(phaseResults(thisSubject).teleporter(thisTele).depths)
            for thisAnalysis = 1:length(phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results)
                for thisFreq = 1:length(phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).GoodFreq)
                    
                    frequency  = phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).GoodFreq(thisFreq);
                    numContigP = round(rayleighCycles / frequency * samplingRate);
                    
                    for thisChan = 1:size(phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).SignifP, 1)
                        pData = squeeze(phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).SignifP(thisChan, thisFreq, :));
                        signifIntervals = identifySignifIntervals(pData, numContigP);
                        if isempty(signifIntervals) == 0
                            for i = 1:size(signifIntervals, 1)
                                outputCell(end + 1, :) = {sessionInfo(thisSubject).subjectID, ...
                                    phaseResults(thisSubject).teleporter(thisTele).name, ...
                                    phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).chanList(thisChan), ...
                                    phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).AnalysisName, ...
                                    phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).GoodFreq(thisFreq), ...
                                    phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).TimeMs(signifIntervals(i, 1)), ...
                                    phaseResults(thisSubject).teleporter(thisTele).depths(thisDepth).results(thisAnalysis).TimeMs(signifIntervals(i, 2))};
                            end
                        end
                    end
                end
            end
        end
    end
end

cell2csv(savePath, outputCell)