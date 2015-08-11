%% make info EVs for the beads task. 

% basic script copied from makeTrialEVs... which was designed to get trial
% onset EVs for the helicopter task. 



codeDir='/Users/mattnassar/matt_work_stuff/Matt/m_files/beadsTask/analysisScripts_MRN';
behavDir='~/Dropbox/BeadsTaskCode/'

addpath(genpath(codeDir), genpath(behavDir));
cd(codeDir);



path(path, '/Users/mattnassar/matt_work_stuff/Matt/workSpaces/heliBehavData')
load data_assembled.mat
outcomeDuration=1;

% for actual fMRI regressors
%load timingData  OLD 16 sub data
ids={d.id}  % select example subject
ind=strmatch('bb914', ids) 
allScanTimes=1.25:2.5:2.5.*226

blockTRs=[ones(length(allScanTimes), 1); ...
    ones(length(allScanTimes), 1).*2; ...
    ones(length(allScanTimes), 1).*3; ...
    ones(length(allScanTimes), 1).*4;]; 

blockTrials=[1:120; 121:240; 241:360; 361:480];
allZero=zeros(length(allScanTimes), length(blockTrials));
trialEvs=struct;

for k = 1:length(d)
    % select subject
    tDat=d(k).mrMain.outcomeTimes;
    id=d(k).id;
    % this is a hack, but not all subjects have 120 trials per block...
    tNum=[d(k).mrMain.statusData.blockCompletedTrials];
    begBlock=find(tNum==1);
    endBlock= [begBlock(2:end)-1 length(tNum)];
    allRegs=nan(length(allScanTimes).*4, length(tNum));
    for j = 1:size(blockTrials, 1)
        % select block
        blockSel=begBlock(j):endBlock(j); % select example block

        %% convolve times with gamma function
        bDat=tDat(blockSel);
        effectSpread=zeros(length(allRegs), length(bDat));
        for f = 1:length(bDat)
            effectTime=bDat(f);
            effectSpread(blockTRs==j,f)=fast_fslgamma(allScanTimes-effectTime,6,3);
        end
        % z-score for final regression coefficients
        allRegs(:,blockSel)=zscore(effectSpread);
    end
     
    eval(sprintf('trialEvs.%s=allRegs;', id));
    
end

save trialEvs.mat trialEvs

