% Call smallMatt.m, named in response to smallArthur.m, a mysteriously
% named piece of code that computes variance inflation in regression
% coeficients for a simulated regression model of BOLD data.



cd /Users/mattnassar/Dropbox/BeadsTaskCode/beadsTask



% 1) get some trials
trialData=getSimTrials(1)
numC=length(trialData.hiLoDiffs)
% 2) set some timing variables:


infoChoiceTime      = 2;
postChoiceJitter    = 3;
firstInfoTime       =.2;
postInfoJitter      = 3;
sideChoiceTime      = 2;
feedbackTime        = 1;
endTrialJitter      = 3;
maxDrawTime         = 5;


totDetTime=infoChoiceTime+firstInfoTime+sideChoiceTime+feedbackTime+maxDrawTime;

% make lists of jitter values
allPCJitter= 2.*postChoiceJitter.*rand(numC,1);
allFIJitter= 2.*postInfoJitter.*rand(numC,1);
allETJitter= 2.*endTrialJitter.*rand(numC,1);
totTime=allPCJitter+allFIJitter+allETJitter;

% these variables will actually depend on subject, but we can model some
% stupid subject behavior just for kicks.
infRespTime =       .5;
drawProb    =       0;
subDrawProb    =       0;



Ps=.1:.1:1;
numReps=20;
clear numSubjs uniqueVar
for i =1:length(Ps)    
    P                   = Ps(i)
    
    for r = 1:numReps
    
    [numSubjs(:,i,r), ~, uniqueVar(:,i, r), names] = smallMatt(trialData,totTime,totDetTime,infoChoiceTime, ...
        firstInfoTime,sideChoiceTime,feedbackTime,maxDrawTime, ...
        allPCJitter,allFIJitter,allETJitter,infRespTime, P);
    end
    
end


% % plot VIF and unique variance for info timing regressor.
% plot(Ps, nanmedian(numSubjs(2,:, :), 3), '.')


defaultPlotParameters
plot(Ps, nanmedian(uniqueVar(2,:, :), 3), '.')
ylabel('unique variance')
xlabel('draw fraction')
title('info timing regressor')
saveas(gcf, 'infoTimingReg_uniqueVar.eps', 'epsc2')
close all


% plot VIF and unique variance for kl divergence.
% plot(Ps, nanmedian(numSubjs(9,:, :), 3), '.')
plot(Ps, nanmedian(uniqueVar(9,:, :), 3), '.')
ylabel('unique variance')
xlabel('draw fraction')
title('info timing regressor')
saveas(gcf, 'klDiverg_uniqueVar.eps', 'epsc2')
close all



% plot VIF and unique variance for kl divergence.
plot(Ps, nanmean(numSubjs(8,:, :), 3), '.')
plot(Ps, nanmean(uniqueVar(8,:, :), 3), '.')





