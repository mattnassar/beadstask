function [trialList]=getSimTrials(blockNum)


blockVars.maxUrnProb  = [.6];
blockVars.hiValSide   = [0 1];

blockList=makeAllCombos(blockVars);
blockOrder=randperm(length(blockList.hiValSide));
blockList=straightStruct(blockList);
changingVars=fieldnames(blockList);

numReps=1;
varList=struct;
varList.hiLoDiffs= [60];
varList.avCorrAvErrDiffs=[40];
varList.infCost=[.01];
varList.expVals=[270 70];
varList.hiValSide=[1];
varList.maxUrnProb=[.6];
varList.startDiff=[-10:2:-2 0 2:2:10];
varList.maxEvenBeads=[20 60];
varList.scaleFactor=[1];

% insert block manipulated variables:
for j=1:length(changingVars)
    eval(sprintf('varList.%s = blockList.%s(blockNum)', changingVars{j}, changingVars{j}));
end

trialList=makeInBeadTrialList(varList, numReps);
