function [numSubjs, uniqueVar, regData, allRegModVars] = getRunEfficiency(trialData, timingData)

%% compute the theoretical efficiency of a run in terms of unique variance in the BOLD
%  response that would result from a given session. 

infoData=getInfoValueFromBlock(trialData)


%  WILL NEED THIS STUFF:  (trialData,totTime,totDetTime,infoChoiceTime,firstInfoTime,sideChoiceTime,feedbackTime,maxDrawTime,allPCJitter,allFIJitter,allETJitter,infRespTime)

%modTrialInfValue; % info value
%modExpMutInfGain; % info quantity
%modTrialStateValue; %overall value of current state
%modDrawActualMutInfoGain; % how much info is actually received at feedback
%modDrawActualValueGain; % how much actual value is realized?
%modDrawActualKLdivergence;

choiceNames={'modTrialInfValue', 'modExpMutInfGain', 'modTrialStateValue'};
infoNames={'modDrawActualMutInfoGain', 'modDrawActualValueGain', 'modDrawActualKLdivergence'};   % 'house'

choicePhaseVars_raw=[infoData.modTrialInfValue infoData.modExpMutInfGain infoData.modTrialStateValue];
choicePhaseVars=zscore(choicePhaseVars_raw);
infoPhaseVars_raw=[infoData.modDrawActualMutInfoGain infoData.modDrawActualValueGain infoData.modDrawActualKLdivergence];   %  infoData.infoType


% mean center choice phase modulator regressors:
choicePhaseVars=choicePhaseVars-repmat(nanmean(choicePhaseVars), length(choicePhaseVars), 1); 
% info phase regressors can only be mean centered after we determine which
% info phases will be observed.

    
ll=length(infoData.modTrialInfValue)
t=0; 

allTrialDraws=nan(ll, 1);
infoOn=timingData.infoOn; 
choiceOn=timingData.choiceOn; 
decTimeOn=timingData.betOn; 
infoOff=timingData.infoOff;
choiceOff=timingData.choiceOff;
decTimeOff=timingData.betChoiceTime;
%expDuration=timingData.feedbackOff-timingData.trialStart;


%% mean center info modulator values JUST on info trials
infoTrials=isfinite(infoOn);
infoPhaseVars=nan(size(infoPhaseVars_raw));
infoPhaseVars(infoTrials,:)=(infoPhaseVars_raw(infoTrials,:)-repmat(nanmean(infoPhaseVars_raw(infoTrials,:)), sum(infoTrials), 1) );


% these trial by trial regressors could be plugged into FSL
allRegModVars.info.vals=infoPhaseVars;
allRegModVars.info.names=infoNames;

allRegModVars.choice.vals=choicePhaseVars;
allRegModVars.choice.names=choiceNames;

allRegModVars.bet.vals=trialData.betRight-nanmean(trialData.betRight);
allRegModVars.bet.names={'betRight'};

allRegModVars.info_raw.vals=infoPhaseVars_raw;
allRegModVars.info_raw.names={'info_raw'};

allRegModVars.choice_raw.vals=choicePhaseVars_raw;
allRegModVars.choice_raw.names={'choice_raw'};



t=600
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
infoPhaseVars(~isfinite(infoPhaseVars))=0;


%

for i = 1:(ll)
   
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


   % now do the exact same process, Left Right decision
    for j = decTimeOn(i):inc:decTimeOn(i)+2 %decTimeOff(i)
        ts=allTrs-j;
        newFdbkResponse=fast_fslgamma(ts,6,3);
        trialFeedbackReg=trialFeedbackReg+newFdbkResponse;
    end
end


allNames=[{'ChoiceReg', 'infoReg', 'decTimeOn'}, choiceNames, infoNames];
allConvRegs=[(choiceReg) (infoReg) trialFeedbackReg (choiceModRegs) (feedbackRegs)];
zAllConvRegs=zscore(allConvRegs);

regData.allConvRegs=allConvRegs;
regData.allConvNAmes=allNames;


% 
%  keyboard;
%  
 
 
%keyboard
%Diagnostics.
% look at timing regressors:
figure
hold on
plot(allTrs, (choiceReg), 'g')
plot(allTrs, infoReg, 'r')
plot(allTrs, trialFeedbackReg, 'b')
xlim([200 300])
set(gca, 'box', 'off')

corr([choiceReg, infoReg, trialFeedbackReg])

% 
% 
% legend('trial on', 'info on', 'feedback on')
% 
% 
% saveas(gcf, 'allRegsOverTime.eps', 'epsc2')
% %close all
% % % 
% 
% 
% figure
% rrr=corr([choiceReg(5:end) infoReg(5:end) trialFeedbackReg(5:end)])
% a=imagesc(rrr)
% set(gca, 'clim', [-1 1])
% colorbar
% saveas(gcf, 'corrMatForArthur.eps', 'epsc2')
% %close all


% 
% plot(choiceReg(10:end), trialFeedbackReg(10:end), '.')
% ylabel('FeedbackTimingReg')
% xlabel('choiceTimingReg')
% set(gca, 'box', 'off')
% % saveas(gcf, 'fdbk_choice_corr.eps', 'epsc2')
% close all
% % 
% %keyboard
% [C,LAGS]=xcov(zscore(choiceReg), zscore(trialFeedbackReg), 10, 'unbiased')
% 
% hold on
% plot(LAGS.*tr, C)
% plot([min(LAGS.*tr) max(LAGS.*tr)], [0 0], '--k')
% plot([0 0], [-1 1], '--k')
% ylim([-1 1])
% ylabel('correlation')
% xlabel('time lag (sec)')
% set(gca, 'box', 'off')
% saveas(gcf, 'crossCovariogram.eps', 'epsc2')
% close all
% 



%keyboard

for i = 1: size(allConvRegs, 2)
    ys=false(size(allConvRegs, 2),1);
    ys(i)=true;
    
    [~,~,R,~,STATS] = regress(allConvRegs(:,ys),[ones(size(allConvRegs,1),1) allConvRegs(:,~ys)]);
    uniqueVar(i)=nanstd(R).^2;  % ie unique variance is just the variance in each regressor that can't be explained by the other regressors!
    
    [~,~,R,~,STATS] = regress(zAllConvRegs(:,ys),[ones(size(allConvRegs,1),1) zAllConvRegs(:,~ys)]);
    zUniqueVar(i)=nanstd(R).^2; % while we really want
    % to maximize unique variance, it is useful to look at unique variance with z-scored inputs so all variable are on the same scale
    % but remember: its really the
    % un-normalized term that will matter for
    % our analysis!
    
    VIF(i)=1./(1-STATS(1));  % compute variance inflation factor.
    
end

numSubjs  = (VIF.* 1./(1./11))+1;
end