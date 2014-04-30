%% 1. Set Variables
function [numberofSubjects,allTrialDraws,time] = arthurssimulator(varList,numReps, repetition,infoChoiceTime,postChoiceJitter,firstInfoTime,postInfoJitter,sideChoiceTime,feedbackTime,endTrialJitter,maxDrawTime,infRespTime)
for m = 1:repetition
    trialDataArray{m}=makeInBeadTrialList(varList, numReps);
end

%% 3. calculations and outcomes
var1 = 0;
var2 = 0;
var3 = 0;
var4 = 0;
var5 = 0;
var6 = 0;
var7 = 0;
var8 = 0;
var9 = 0;
tottottime = 0;

for q = 1: length(trialDataArray)
    trialData = trialDataArray{q};
    % make lists of jitter values
    allPCJitter= postChoiceJitter.*rand(44,1);
    allFIJitter= postInfoJitter.*rand(44,1);
    allETJitter= endTrialJitter.*rand(44,1);
    % this is the total duration of the trial before adding jitter
    totDetTime = infoChoiceTime+firstInfoTime+sideChoiceTime+feedbackTime+maxDrawTime;
    totTime=allPCJitter+allFIJitter+allETJitter+totDetTime;
    tottottime = tottottime + sum(totTime)/60;
    [numSubjs,allTrialDraws] = smallMatt(trialData,totTime,totDetTime,infoChoiceTime,firstInfoTime,sideChoiceTime,feedbackTime,maxDrawTime,allPCJitter,allFIJitter,allETJitter,infRespTime);
    var1 = var1+numSubjs(1);
    var2 = var2+numSubjs(2);
    var3 = var3+numSubjs(3);
    var4 = var4+numSubjs(4);
    var5 = var5+numSubjs(5);
    var6 = var6+numSubjs(6);
    var7 = var7+numSubjs(7);
    var8 = var8+numSubjs(8);
    var9 = var9+numSubjs(9);
end
time = tottottime/repetition;
numberofSubjects(1) = var1/repetition;
numberofSubjects(2) = var2/repetition;
numberofSubjects(3) = var3/repetition;
numberofSubjects(4) = var4/repetition;
numberofSubjects(5) = var5/repetition;
numberofSubjects(6) = var6/repetition;
numberofSubjects(7) = var7/repetition;
numberofSubjects(8) = var8/repetition;
numberofSubjects(9) = var9/repetition;
end
