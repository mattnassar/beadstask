%% CODE to check to get unique variance on from a task run



%% get everything we need on the path.
clear classes

rootDir='/Users/mattnassar/matt_work_stuff/Matt/m_files/beadsTask/'
cd(rootDir);
addpath(genpath(rootDir))



%% load subject data
subFileDir='/Users/mattnassar/Dropbox/BeadsTaskCode/BeadsTask5data/';

% DATA FROM OUT OF SCANNER PILOT FOR FORCED DRAW VERSION:
%subNames={'BE1829',	'KK1101',	'testbug'};
subNames={'TQ3543'};


% DATA FROM INITIAL fMRI PILOT:
% subNames={'DZ572',	'HE3515',	'KQ3513',	'LT3516',	'LZ3512',	'NN3511', ...
%   'OK547',	'SI3517',   'EO3518', 	'HU2138',	'LL2026',	'LX1803', ...
%   'MC1906',	'NS969',	'QB2027',	'TA3514'};

% loop through each subject and load data, store in "allSubjData" structure
% and save that structure so we don't need to use this slow loading code
% again.
for i = 1:length(subNames)
    fileName=fullfile(subFileDir, subNames{i});
    [allData]=beadTaskDataLoader5(fileName);
    
    % put data structures in a bigger structure including all subjects.
    eval(sprintf('allSubjData.%s=allData;', subNames{i}));
    
    % beadTaskDataLoader changes directories, so lets go back manually.
    cd(subFileDir)
end


%% Get an example block and check the block efficiency. 


[dataStruct]=unpackBeadsData(allSubjData.TQ3543.realScan2.block2.statusData)
[numSubjs, uniqueVar, timingData, trialData, regData, allRegModVars]=checkBlockEfficiency(dataStruct)

% there seems to be a consistent negative correlation between "info" and
% "choice" phases... also a less prominent negative correlation between "info" and
% "feedback" phaes.







