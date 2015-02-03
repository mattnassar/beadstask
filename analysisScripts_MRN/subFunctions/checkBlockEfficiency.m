function [numSubjs, uniqueVar, timingData, trialData, regData, allRegModVars]=checkBlockEfficiency(dataStruct)

% This function serves a few purposes. 

% 1) it takes a raw dataStruct and extracts timing data and puts it in a
% timingData structure.

% 2) extracts trialData and puts it in a trialData structure

% 3) it uses both of these things to compute model efficiency. 





ll=length(dataStruct.infoButtonSide);
trialData.maxUrnProb=dataStruct.curr_maxUrnProb;
trialData.startLeft=dataStruct.curr_start_tokensLeft;
trialData.startRight=dataStruct.curr_start_tokensRight;
trialData.infCost=dataStruct.curr_infCost;
trialData.hiValue=max([dataStruct.curr_rewCorrLeft dataStruct.curr_rewCorrRight], [], 2);
trialData.loValue=min([dataStruct.curr_rewCorrLeft dataStruct.curr_rewCorrRight], [], 2);
trialData.inValue=dataStruct.curr_penErrLeft;
trialData.curr_trialInfList=dataStruct.curr_trialInfList(1:ll);
trialData.draw=isfinite(dataStruct.infoOn);
trialData.hiValSide=dataStruct.curr_rewCorrRight>dataStruct.curr_rewCorrLeft;
trialData.betRight=dataStruct.curr_choice==2;
trialData.infoButtonSide=dataStruct.infoButtonSide;

%keyboard
inT=dataStruct.startTime;
%keyboard
timingData=struct;
timingData.infoOff=dataStruct.infoOn-inT;      % this was event was mistakenly named
% if isfield(dataStruct, 'realInfoOn')
%     timingData.infoOn=dataStruct.realInfoOn-inT;
% else
    % if we don't have a timestamp, compute it.
    disp('Computing infoOn timestamp based on info off')
    timingData.infoOn  =  timingData.infoOff-1.0776;
% end
timingData.choiceOn      =dataStruct.preChoiceOn-inT;
timingData.feedbackOn    =dataStruct.fdbkOn-inT;     % now defunct bc there is no feedback
timingData.choiceOff     =dataStruct.preChoiceOff-inT;
timingData.feedbackOff   =dataStruct.fdbkOff-inT;   % this is now defunct bc there is no feedback
timingData.trialStart    =dataStruct.trialStartTime-inT;
timingData.betOn         =dataStruct.betOn-inT;
timingData.tooSlow       =dataStruct.tooSlowOnTime;

% new timestamps added by Arthur ~ 07/14/14
timingData.betOn2        =dataStruct.enterDecTime-inT; % ok, this occurs right before betOn
timingData.freeDrawOn    =dataStruct.freeDrawOn-inT;
timingData.freeDrawOff   =dataStruct.freeDrawOff-inT; % this is the same as betOn
timingData.betChoiceTime =dataStruct.betchoicetime-inT; % when subject actually presses button indicating bet



if isfield(dataStruct, 'infChoiceTime')
    timingData.infoChoiceMade=dataStruct.infChoiceTime-inT;
else
    timingData.infoChoiceMade=nan(size(timingData.choiceOn));
end


for i = 1:ll
    
    if isfield(dataStruct, 'extraDrawTime')
        timingData.extraButtonPush{i}=dataStruct.extraDrawTime{i}-unique(inT);
    else
        timingData.extraButtonPush{i}=[];
    end
end







%getInfoValueFromModel
[numSubjs, uniqueVar, regData, allRegModVars] = ...
    getRunEfficiency(trialData, timingData);


