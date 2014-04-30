function [numSubjs, allTrialDraws] = smallMatt(trialData,totTime,totDetTime,infoChoiceTime,firstInfoTime,sideChoiceTime,feedbackTime,maxDrawTime,allPCJitter,allFIJitter,allETJitter,infRespTime)
getInfoValueFromMode;
%modTrialInfValue; % info value
%modExpMutInfGain; % info quantity
%modTrialStateValue; %overall value of current state
%modDrawActualMutInfoGain; % how much info is actually received at feedback
%modDrawActualValueGain; % how much actual value is realized?
%modDrawActualKLdivergence;

choiceNames={'modTrialInfValue' 'modExpMutInfGain' 'modTrialStateValue'};
infoNames={'modDrawActualMutInfoGain' 'modDrawActualValueGain' 'modDrawActualKLdivergence'};

choicePhaseVars=zscore([modTrialInfValue modExpMutInfGain modTrialStateValue]);
infoPhaseVars=[modDrawActualMutInfoGain modDrawActualValueGain modDrawActualKLdivergence];

% mean center choice phase modulator regressors:
choicePhaseVars=choicePhaseVars-repmat(nanmean(choicePhaseVars), length(modTrialInfValue), 1); 
% info phase regressors can only be mean centered after we determine which
% info phases will be observed.


% run through a bunch of "subject types" ie drawers and non-drawers

% compute bead difference for each trial.
if ~unique(trialData.hiValSide)
    beadDifference=trialData.startLeft-trialData.startRight;
else
    disp('put high value on left')
end
    
    
t=0; allTrialDraws=nan(length(modTrialInfValue), 1);
infoOn=nan(length(modTrialInfValue), 1); infoOff=nan(length(modTrialInfValue), 1);
choiceOn=nan(length(modTrialInfValue), 1); choiceOff=nan(length(modTrialInfValue), 1);
feedbackOn=nan(length(modTrialInfValue), 1); feedbackOff=nan(length(modTrialInfValue), 1);

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
    if todrawornot(beadDifference(i));
        t=t+allPCJitter(i);
        infoOn(i)=t;
        t=t+firstInfoTime;
        infoOff(i)=t;
        t=t+allFIJitter(i);
        trialDraws=trialDraws+1;
        inT=t;
        
        % reset the bead difference according to initial draw
        if trialData.corrUrnSide(i)==1;
            beadDifference(i)=beadDifference(i) +  2.*((rand<.6)-.5);
        else
            beadDifference(i)=beadDifference(i) +  2.*((rand<.4)-.5);
        end

        % check for second draw... 
        while todrawornot(beadDifference(i)) && (t-inT)<maxDrawTime
            % if we draw another bead, then reset bead difference
            if trialData.corrUrnSide(i)==1;
                beadDifference(i)=beadDifference(i) +  2.*((rand<.6)-.5);
            else
                beadDifference(i)=beadDifference(i) +  2.*((rand<.4)-.5);
            end   
            % and waste some time
            trialDraws=trialDraws+1;
            t=t+infRespTime;
        end
    else
        t=t+allETJitter(i);
    end
    
    t=t+sideChoiceTime;
    
    feedbackOn(i)=t;
    t=t+feedbackTime;
    feedbackOff(i)=t;
    
    t=t+allETJitter(i);
    allTrialDraws(i)=trialDraws;
    totTrialT(i)=t-inTrialT;
    totTimeSaved(i)=expDuration(i)-totTrialT(i);
end
extraTrials=floor(sum(totTimeSaved)./nanmean(expDuration));


%% mean center info modulator values JUST on info trials
infoTrials=isfinite(infoOn);
infoPhaseVars(infoTrials,:)=(infoPhaseVars(infoTrials,:)-repmat(nanmean(infoPhaseVars(infoTrials,:)), sum(infoTrials), 1) );


tr=2.5;  % how long will it take to scan the entire brain once
numTrs=ceil(t./tr)+5 ; % how many TRs will we need total to cover the entire task block
allTrs=tr/2:tr:numTrs.*tr;  %lets list them.
inc=.5;

% ok, lets loop through trials
choiceReg=zeros(size(allTrs))';
infoReg=zeros(size(allTrs))';
trialFeedbackReg=zeros(size(allTrs))';
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




   % now do the exact same process, Feedback
    for j = feedbackOn(i):inc:feedbackOff(i)
        ts=allTrs-j;
        newFdbkResponse=fast_fslgamma(ts,6,3);
        trialFeedbackReg=trialFeedbackReg+newFdbkResponse;
    end
end

allNames=[{'ChoiceReg', 'infoReg', 'feedbackOn'}, choiceNames, infoNames];
allConvRegs=[(choiceReg) (infoReg) trialFeedbackReg (choiceModRegs) (feedbackRegs)];
zAllConvRegs=zscore(allConvRegs);

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

numSubjs  = (VIF.* 1./(1./11))+1;
end