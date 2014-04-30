%% Beads task scanner version mock-up...

% the goal of this script is to create a mock version of the beads task
% that could be run in the scanner. This will involve creating a bunch of
% trial types, blocks, etc.

% then we will want to examine what the BOLD signals related to our
% variables of interest will look like (ie how much unique variance related
% to the expected value of information would be in our regression model?
% Could the amount of variance be increased by adding jitter/delays between
% trials? how well can we dissociate value from quantity of information?)


% steps
% 1) get some trials
% 2) get regression variables per trial
%   A) value of information (model, or side, or subjective)
%   B) quantity of information (expected number of bits)
%   C) expected value of state (from model)
%   D) feedback variables: mut inf gain, delta V
% 3) assign timestamps to trials
% 4) model bold response for trials according to regression terms
% 5) optimize timing variables



% In order to maximize our ability to dissociate signals related to
% each variable of interest, you can change several things:

% 1) The specific types of trials that we are choosing.  For example I
% think that increasing the overall value offset in the expected value
% manipulation condition should decrease the negative correlations between 
% expected value and expected "value of information", and thus increase the
% unique variance related to each term.

% 3) Yeah, i skipped 2. I'm following the numbering scheme above. So we're
% talking about timestamps. This is important because the BOLD response is
% quite slow, so unfortunately we wont get an exact measure of the trial
% variables listed above.  Instead, we'll get a smeared mix of each
% measurement along with the things that happened before and after that
% measurement. The goal here is to choose the timestamps that maximize the
% amount of unique variance that we KEEP even after said smearing is
% complete.




%% 1) get some trials
%    This is one thing that you can mess around with to try to optimize our
%    scanner task. As you see trialData gives us a list of trial
%    conditions... if you go into the function getSimTrials you can change
%    the exact set of conditions to see how this effects our overall levels
%    of unique variance (and thus ability to identify neural signals)

%trialData=getSimTrials(1);
clear classes
blockNum = 1;
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
varList.infCost=[.5];
varList.expVals=[140 70];
varList.hiValSide=[1];
varList.maxUrnProb=[.6];
varList.startDiff=[-10:2:10];
varList.maxEvenBeads=[20 30];
varList.scaleFactor=[1];

% insert block manipulated variables:
for j=1:length(changingVars)
    eval(sprintf('varList.%s = blockList.%s(blockNum)', changingVars{j}, changingVars{j}));
end

trialData=makeInBeadTrialList(varList, numReps);

%% 2) get variables

% A) value of information

getInfoValueFromMode
modTrialInfValue; % info value
modExpMutInfGain % info quantity
modTrialStateValue %overall value of current state
modDrawActualMutInfoGain % how much info is actually received at feedback
modDrawActualValueGain % how much actual value is realized?
modDrawActualKLdivergence



%  One factor limiting our power will be correlations between our variables
%  of interest across all trials. You can look at this before even thinking 
%  about timing in order to tell how well dissociated our terms are by
%  trial

%   The code also spits out a heatmap of correlations... diagonal should be
%   red indicating high correlation, but other correlations should be
%   fairly moderate (ie absolutely small).

rho=corr([modTrialInfValue modExpMutInfGain modTrialStateValue]);
surf2d(rho)
colorbar
set(gca, 'clim', [-1 1])
saveas(gcf, 'rawCorrMatrix.eps', 'epsc2')
close all

%  This task also dissociates some important factors during the information
%  presentation part of the task. Looking at correlations can give you a
%  sense of how good (or not good) a job it does at this.

%  However, there is a bit of a problem with this approach: The subject
%  will not see information on each trial. So a better way to do this would
%  be to use our pilot data to simulate some realistic behavior and
%  determine how well we do when we have limited draw trials. It will
%  certainly be worse, but i'm not sure whether it will be catostrophic.

choiceNames={'modTrialInfValue' 'modExpMutInfGain' 'modTrialStateValue'};
infoNames={'modDrawActualMutInfoGain' 'modDrawActualValueGain' 'modDrawActualKLdivergence'};

choicePhaseVars=zscore([modTrialInfValue modExpMutInfGain modTrialStateValue]);
infoPhaseVars=[modDrawActualMutInfoGain modDrawActualValueGain modDrawActualKLdivergence];

rho=corr(infoPhaseVars);
surf2d(rho)
set(gca, 'clim', [-1 1])
colorbar
close all


%  OK, this chunk of code makes a plot that is on my poster... the point is just
%  to demonstrate HOW the decision phase variables are decorrelated.
%  Should be pretty intuitive once you see it...


hold on
modTrialInfValue=modTrialInfValue.*trialData.scaleFactor;
getCbColors
plot([0 0], [0 1], '--k')
startDiff=trialData.startRight-trialData.startLeft;
[xes, i]=sort(startDiff);

a=plot(-xes, unitNorm(modTrialInfValue(i)), '-', 'color', cbColors(4,:));
b=plot(-xes, unitNorm(modExpMutInfGain(i)), '-r', 'color', cbColors(3,:));

highValTrials=trialData.hiValue==max(trialData.hiValue);
stVal=modTrialStateValue(highValTrials);
[x2, i]=sort(startDiff(highValTrials));
c=plot(-x2, (stVal(i)-min(modTrialStateValue))./(max(modTrialStateValue)-min(modTrialStateValue)), 'color', cbColors(5,:));

lowValTrials=trialData.hiValue==min(trialData.hiValue);
stVal=modTrialStateValue(lowValTrials);
[x2, i]=sort(startDiff(lowValTrials));
y=unitNorm(stVal);
plot(-x2, (stVal(i)-min(modTrialStateValue))./(max(modTrialStateValue)-min(modTrialStateValue)), 'color', cbColors(5,:))

ylabel('Normalized values')
xlabel('Bead difference (high-low)')
aa=legend([a b c], 'Info value', 'Info quantity', 'State value');
set(aa, 'box', 'off', 'location', 'east')
set(gca, 'box', 'off')
saveas(gcf, 'preConvolvedData.eps', 'epsc2')
close all


% not sure what i was doing here...
%plot(modExpMutInfGain./max(modExpMutInfGain), modTrialInfValue./max(modTrialInfValue), '.')



%% 3) assign timestamps to trials

% ok, the goal here is to simulate an entire session across time. In order
% to do this we need to know how long each event will take. Actually, it would be 
% good to know what the current settings are for a bunch of these
% variables...


% trialsetup, all in seconds


infoChoiceTime      = 3;
postChoiceJitter    = 8;
firstInfoTime       =.2;
postInfoJitter      = 0;
sideChoiceTime      = 2;
feedbackTime        = 2;
endTrialJitter      = 8;
maxDrawTime         = 2;

% this is the total duration of the trial before adding jitter
totDetTime=infoChoiceTime+firstInfoTime+sideChoiceTime+feedbackTime+maxDrawTime;


% make lists of jitter values
allPCJitter= postChoiceJitter.*rand(length(modTrialInfValue),1);
allFIJitter= postInfoJitter.*rand(length(modTrialInfValue),1);
allETJitter= endTrialJitter.*rand(length(modTrialInfValue),1);
totTime=allPCJitter+allFIJitter+allETJitter+totDetTime;

% these variables will actually depend on subject, but we can model some
% stupid subject behavior just for kicks.

% actually, the above comment was from when i first wrote this code. But
% now i think it makes sense to plug in our pilot data from the fMRI
% version.


infRespTime =       .5;   % how long does it take the subject to respond
drawProb    =       0;     % these two terms are currently over-ridden by the 
subDrawProb    =       0;   % loop below. Right now we are simulating a fixed  
                              % probability of drawing on each trial... but
                               % in reality we should be computing it
                               % separately according to the conditions and
                               % according to our pilot data. 

drawProbs=0:.1:1;

% run through a bunch of "subject types" ie drawers and non-drawers
for k = 1:length(drawProbs)

    drawProb=drawProbs(k); % probability of drawing on the first go
    subDrawProb=drawProbs(k);    % probability of drawing subsequent bead
    t=0; allTrialDraws=nan(length(modTrialInfValue), 1);
    infoOn=nan(length(modTrialInfValue), 1);
    choiceOn=nan(length(modTrialInfValue), 1);
    % run through each trial and compute the timing of some important
    % events
    for i = 1:length(modTrialInfValue)
        inTrialT=t;
        expDuration(i)=totTime(i)+totDetTime;
        choiceOn(i)=t;
        t=t+infoChoiceTime;
        choiceOff(i)=t;
        
        % also simulate some behavior. 
        % ie. will subject draw?
        
        
        trialDraws=0;
        if rand<drawProb
            t=t+allPCJitter(i);
            infoOn(i)=t;
            t=t+firstInfoTime;
            infoOff(i)=t;
            t=t+allFIJitter(i);
            trialDraws=trialDraws+1;
            x=rand; inT=t;
            while x<subDrawProb & (t-inT)<maxDrawTime
                trialDraws=trialDraws+1;
                t=t+infRespTime;
                x=rand;
            end
        else
            t=t+allETJitter(i);
        end
        
        t=t+sideChoiceTime;
        t=t+feedbackTime;
        t=t+allETJitter(i);
        allTrialDraws(i)=trialDraws;
        totTrialT(i)=t-inTrialT;
        totTimeSaved(i)=expDuration(i)-totTrialT(i);
    end
    extraTrials(k)=floor(sum(totTimeSaved)./nanmean(expDuration));
    
end



% how many extra trials could we fit for each level of drawing?
plot(drawProbs, extraTrials)
ylabel('extra trials')
xlabel('draw probability')
saveas(gcf, 'trialsVsDrawProb.eps', 'epsc2')
close all



% OK, now we have 1) task event timings, 2) event modulator values... so we
% we know what signals we'd expect to see at the level of neural activity.
% unfortunately we're measuring the BOLD response, which is slow and spread
% in time. In order to know what our signals would look like at the level
% of the BOLD response we need to convolve our expected neural signal with 
% the haemodynamic response function HRF (ie the expected BOLD response to a
% brief neural activation). 

% one function that is used to approximate the HRF is a gamma function
% fast_fslgamma makes gamma convolved with stick. need to convolve with
% box in order to model choice and feedback epochs



% t./60   % not sure what i was doing here..
tr=2.5;  % how long will it take to scan the entire brain once
numTrs=ceil(t./tr)+5 ; % how many TRs will we need total to cover the entire task block
allTrs=tr/2:tr:numTrs.*tr;  %lets list them.
inc=.5;

% ok, lets loop through trials
choiceReg=zeros(size(allTrs))';
infoReg=zeros(size(allTrs))';
feedbackRegs=zeros(length(allTrs), size(infoPhaseVars, 2));
choiceModRegs=zeros(length(allTrs), size(choicePhaseVars, 2));
for i = 1:length(modTrialInfValue)
    
    % run through choice phase (by steps of inc) and simulate BOLD responses
    trialReg=zeros(size(allTrs))';
    for j = choiceOn(i):inc:choiceOff(i)
        ts=allTrs-j;
        newResponse = fast_fslgamma(ts,6,3);
        choiceReg=choiceReg+newResponse;
        trialReg=trialReg+newResponse;  % regressor for this phase of the trial
    end
    
    trialReg=trialReg./length(choiceOn(i):inc:choiceOff(i));
    
    %    plot(repmat(choicePhaseVars(i,:), length(allTrs), 1).*repmat(trialReg, 1, size(choicePhaseVars, 2)))
    
    
    choiceModRegs=choiceModRegs+repmat(choicePhaseVars(i,:), length(allTrs), 1).* ...
        repmat(trialReg, 1, size(choicePhaseVars, 2));   % Modulator terms (ie Trial time regressor * trial variable regressor)
    
    
    % now do the exact same process, but run through the info phase
    trialReg=zeros(size(allTrs))';
    for j = infoOn(i):inc:infoOff(i)
        ts=allTrs-j;
        newInfoResponse=fast_fslgamma(ts,6,3);
        infoReg=infoReg+newInfoResponse;
        trialReg=trialReg+newInfoResponse; % compute an info phase BOLD response
        
    end
    trialReg=trialReg./length(choiceOn(i):inc:choiceOff(i));
    feedbackRegs=feedbackRegs+repmat(infoPhaseVars(i,:), length(allTrs), 1).* ...
        repmat(trialReg, 1, size(infoPhaseVars, 2)); % multiply it by trial variables to get info phase modulator resopnses
    
end



% look at expected responses of areas responeding to choice phase (red) or
% info phase (green). 
hold on
plot(zscore(choiceReg), 'r')
plot(zscore(infoReg), 'g')
xlim([100 200])
ylabel('phase regressors')
xlabel('time')
saveas(gcf, 'demoRegs.eps', 'epsc2')
close all



% in theory you could plot everything out this way in time, but the rest of
% the regressors will look a lot less intuitive.
% so instead lets just look at the correlations between different
% regressors, this will give you a sense of which factors might be
% correlated in the current situation. the allNames variable will tell you
% which variables are being plotted in the correlation matrix.


allNames=[{'ChoiceReg', 'infoReg'}, choiceNames, infoNames];
allConvRegs=[(choiceReg) (infoReg) (choiceModRegs) (feedbackRegs)];
zAllConvRegs=zscore(allConvRegs);
rho=corr(zAllConvRegs);
surf2d(rho)
colorbar
set(gca, 'clim', [-1 1])
saveas(gcf, 'allRegCorrMatrix.eps', 'epsc2')
close all



% This correlation business might seem a bit weird... our endgame is that
% we want to be able to uniquely identify signals that are specific to one
% thing, not signals that are shared. Correlations hurt us, but it would be
% nicer to have a direct measure of how much variance we DO have available
% to us. So here we'll do it by computing the unique variance for each term
% and the variance inflation factor which quantifies the extent to which
% our estimation variance will be amplified due to input correlations. 

% calculate unique variance for each regressor
for i = 1: size(allConvRegs, 2)
    ys=false(size(allConvRegs, 2),1);
    ys(i)=true;
    
    [B,BINT,R,RINT,STATS] = regress(allConvRegs(:,ys),[ones(size(allConvRegs,1),1) allConvRegs(:,~ys)]);
    uniqueVar(i)=nanstd(R).^2;  % ie unique variance is just the variance in each regressor that can't be explained by the other regressors!
    
    [B,BINT,R,RINT,STATS] = regress(zAllConvRegs(:,ys),[ones(size(allConvRegs,1),1) zAllConvRegs(:,~ys)]);
    zUniqueVar(i)=nanstd(R).^2; % while we really want 
                                % to maximize unique variance, it is useful to look at unique variance with z-scored inputs so all variable are on the same scale
                                % but remember: its really the
                                % un-normalized term that will matter for
                                % our analysis!
                               
    VIF(i)=1./(1-STATS(1));  % compute variance inflation factor.
    
end


%% make some plots:

% the key variable is how much unique variance we have in each term. 
bar(1:size(allConvRegs, 2), uniqueVar, .25)
ylabel('unique variance')
xlabel('regressor number')

%However, its a bit tough to look at those, bc they are un-normalized.
%another thing to look at is how much your variance will be inflated due to
%correlations
bar(1:size(allConvRegs, 2), VIF, .25)
ylabel('variance inflation')
xlabel('regressor number')
saveas(gcf, 'VIF.eps', 'epsc2')


%
%     %% lets calculate the number of subjects that would be required to
%     %% produce the same estimate variance as 12 subjects with orthogonal
%     %% regressors:
%
numSubjs  = (VIF.* 1./(1./11))+1
sum(totTime)/60
bar(1:size(allConvRegs, 2), numSubjs, .25)
ylabel('N equivalent to 12')
xlabel('regressor number')


    
