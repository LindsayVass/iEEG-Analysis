% Create a GUI that asks the user for information about this patient

% Fields
clear P;
P.PatientNumber = {'' , 'Patient Number'};
P.Age = { '' [1 100] };
P.Gender = { 'M|F' };
P.Race = { { '{}', 'American Indian/Alaska Native', 'Asian', 'Black/African American', 'Native Hawaiian/Pacific Islander', 'White/Caucasian'} };
P.Ethnicity = {'Hispanic|Not Hispanic'};
P.BehavioralIssues = { 'No reported behavioral or memory issues', 'Behavioral Issues' };
P.ImplantDate = { [date], 'Implant Date' };
P.ElectrodeType = { 'sEEG', 'Electrode Type' };
P.ElectrodeCoverage.LHipp = { {'{0}' '1'} };
P.ElectrodeCoverage.RHipp = { {'{0}' '1'} };
P.ElectrodeCoverage.LFrontal = { {'{0}' '1'} };
P.ElectrodeCoverage.RFrontal = { {'{0}' '1'} };
P.ElectrodeCoverage.LParietal = { {'{0}' '1'} };
P.ElectrodeCoverage.RParietal = { {'{0}' '1'} };
P.ElectrodeCoverage.LTemporal = { {'{0}' '1'} };
P.ElectrodeCoverage.RTemporal = { {'{0}' '1'} };
P.ElectrodeCoverage.LOccipital = { {'{0}' '1'} };
P.ElectrodeCoverage.ROccipital = { {'{0}' '1'} };
P.ElectrodeCoverage.Other = { '' };
P.ImplantHistory = { 'No previous implant', 'Implant History' };
P.ResectionHistory = { 'No previous resection', 'Resection History' };

Patient = StructDlg(P, 'Patient Information');
