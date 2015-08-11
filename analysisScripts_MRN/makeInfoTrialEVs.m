function trialEVs=makeInfoTrialEVs(trialTimes, TR, numTR, blkNum, subName, outFolder)

% this function should take scan sequence parameters (TR, numTR) and
% event timings and duration to create a set of mumford and poldrack
% style trial regressors.

% we'll model trial as an instantaneous event (this is how i did it for the 
% the helicopter task, so i'm just going to try to keep things consistent. 

% NOTE: one difference between this script and the helicopter one is that
% the heli trialEV script z-scored all regressors. This doesn't seem
% exactly right... so here I'm going to leave them un-normalized and do
% normalization (z-scoring) after summing over the non-relevant trials. 



numBlocks=max(blkNum);
allScanTimes=TR./2:TR:numTR.*TR;
blockTRs=[ones(length(allScanTimes), 1); ...
    ones(length(allScanTimes), 1).*2; ...
    ones(length(allScanTimes), 1).*3; ...
    ones(length(allScanTimes), 1).*4;]; 

for i =1:numBlocks
    if i == 1
    allSesScanTimes=allScanTimes;
    else
    allSesScanTimes=[allSesScanTimes,    allScanTimes+  (i-1).*numTR.*TR];
    end
end

trialEVs=struct;
allRegs=nan(length(allScanTimes).*4, length(trialTimes));

%% convolve times with gamma function
effectSpread=zeros(length(allSesScanTimes), length(trialTimes));
for f = 1:length(trialTimes)
    effectTime=trialTimes(f);
    effectSpread(:,f)=fast_fslgamma(allSesScanTimes-effectTime,6,3);
    % get rid of "across block
    effectSpread(~(blkNum(f)==blockTRs), f)=0;
end
% z-score for final regression coefficients
allRegs=(effectSpread);

trialEVs.subName=subName;
trialEVs.trialTimes=trialTimes;
trialEVs.TR=TR;
trialEVs.numTR=numTR;
trialEVs.blkNum=blkNum;
trialEVs.allRegs=allRegs;

outFN=fullfile(outFolder, sprintf('trialInfoEvs_%s.mat',  subName));
save(outFN, 'trialEVs');



