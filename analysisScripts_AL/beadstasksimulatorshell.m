clear classes
clc
blockNum = 1;
blockVars.maxUrnProb  = .6;
blockVars.hiValSide   = [0 1];
blockList=makeAllCombos(blockVars);
blockOrder=randperm(length(blockList.hiValSide));
blockList=straightStruct(blockList);
changingVars=fieldnames(blockList);
numReps=1;
varList=struct;
varList.hiValSide=1;
varList.maxUrnProb=.6;
varList.startDiff=-10:2:10;
varList.maxEvenBeads=[20 30];
varList.scaleFactor=1;
varList.expVals = [270,70];
varList.avCorrAvErrDiffs = 40;
varList.hiLoDiffs = 60;
varList.infCost = .5;
infRespTime =.5;
for j=1:length(changingVars)
    eval(sprintf('varList.%s = blockList.%s(blockNum);', changingVars{j}, changingVars{j}));
end

expVals = varList.expVals;
avCorrAvErrDiffs = varList.avCorrAvErrDiffs;
hiLoDiffs = varList.hiLoDiffs;
infCost = varList.infCost;

%%

counter = 0;
trialList = [];
for infoChoiceTime = 2:1:2
    for postChoiceJitter = 6:1:6
        for firstInfoTime = .2:.2:.2
            for postInfoJitter = 6:1:6
                for sideChoiceTime = 2:1:2;
                    for feedbackTime = 1:1:1;
                        for endTrialJitter = 6:1:6  % 2~9
                            for maxDrawTime = 5:1:5  %3~4
                                trialTime = infoChoiceTime+firstInfoTime+sideChoiceTime+feedbackTime+maxDrawTime+.5*postChoiceJitter + .5*postInfoJitter + .5*endTrialJitter;
                                minutes = (trialTime*44)/60;
                                if minutes <= 15 && minutes >= 14
                                    counter = counter+1;
                                    trialList(counter,:) = [infoChoiceTime, firstInfoTime, sideChoiceTime, feedbackTime, maxDrawTime, postChoiceJitter, postInfoJitter, endTrialJitter];
                                end
                            end
                        end
                    end
                end
            end
            
        end
    end
end

%%
numcondition = size(trialList,1);
allData = struct;
allData.numSubjs = {};
allData.time = [];
allData.timingconditions = {};
allData.allTrialDraws = {};



repetition = 1;

for q = 1: numcondition
    disp(q/numcondition)
    [numberofSubjects,allTrialDraws,time] = arthurssimulator(varList,numReps, repetition,trialList(q,1),trialList(q,6),trialList(q,2),trialList(q,7),trialList(q,3),trialList(q,4),trialList(q,8),trialList(q,5),infRespTime);
    allData.numSubjs{q} = numberofSubjects;
    allData.time(q) = time;
    allData.timingconditions{q} = trialList(q,:);
    allData.allTrialDraws{q} = allTrialDraws;
end
disp('all done')
save allData