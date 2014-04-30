function [allTrs, fslregressor] = getRegressors(statusData)
%% step 1. extract trial data and timing data
allN={statusData.name};
goodStatus=~cellfun(@isempty, allN);
statusData=statusData(goodStatus);
dataStruct.curr_trialInfList={statusData.curr_trialInfList};
allNames=fieldnames(statusData);
for j = 1:length(allNames)
    % get data from a field
    dataStruct.(allNames{j}) = {statusData.(allNames{j})};
end
dataStruct=straightStruct(dataStruct);
for hp = 1:length(dataStruct.drawBegTime)
    if isempty(dataStruct.drawBegTime{hp})
        dataStruct.drawBegTime{hp} = [nan];
    end
end
for j = 1:length(allNames)
    % get data from a field
    lenarray = zeros(length(dataStruct.(allNames{j})),1);
    for q = 1:length(dataStruct.(allNames{j}))
        lenarray(q) = length(dataStruct.(allNames{j}){q});
    end
    uniq = unique(lenarray);
    if sum(uniq) == 1 && ~strcmp(allNames{j}, 'curr_trialInfList')
        dataStruct.(allNames{j}) = cell2mat(dataStruct.(allNames{j}));
    end
end


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
trialData.leftrightbet = dataStruct.curr_choice;

t=600;
tr=.1;  % how long will it take to scan the entire brain once
numTrs=ceil(t./tr)+5 ; % how many TRs will we need total to cover the entire task block
allTrs=0:tr:(numTrs+500).*tr;
inT=dataStruct.startTime(1); %% when the task started by a trigger from the scanner
allTrs = allTrs';
allTrs = [allTrs zeros(length(allTrs),1) zeros(length(allTrs),1) zeros(length(allTrs),1) zeros(length(allTrs),1) zeros(length(allTrs),1) zeros(length(allTrs),1) zeros(length(allTrs),1)];

timingData=struct;
timingData.preChoiceOn = dataStruct.preChoiceOn - inT; %regressor 1 on time
timingData.preChoiceOff = dataStruct.preChoiceOff - inT; %regressor 1 off time/regressor 2 on time

timingData.infChoiceTime = dataStruct.infChoiceTime - inT; %regressor 2 off time
for m = 1:length(timingData.infChoiceTime)
    if ~isfinite(timingData.infChoiceTime(m)) % if the regressor 2 off timestamp is unavailable, it means too slow was turned on
        timingData.infChoiceTime(m) = dataStruct.tooSlowOnTime(m) - inT;
    end
end
timingData.drawBegTime = dataStruct.drawBegTime - inT; %regressor 3 on time
timingData.infoOn=dataStruct.infoOn-inT; %regressor 3 off time

timingData.freeDrawOn = dataStruct.freeDrawOn-inT; %regressor 4 on time
timingData.freeDrawOff = dataStruct.freeDrawOff-inT; %regressor 4 off time

timingData.betOn = dataStruct.betOn - inT; %regressor 5 on time
timingData.fdbkOn = dataStruct.fdbkOn - inT; %regressor 5 off time/ regressor 6 on time
timingData.fdbkOff = dataStruct.fdbkOff - inT; %regressor 6 off time

timingData.tooSlowOnTime = dataStruct.tooSlowOnTime - inT; %regressor 7 on
timingData.tooSlowOffTime = dataStruct.tooSlowOffTime - inT; %regressor 7 off

for m = 1: length(allTrs)
    tmp = allTrs(m,1);
    for n = 1:length(timingData.preChoiceOn)
        if timingData.preChoiceOn(n) <= tmp && tmp < timingData.preChoiceOff(n)
            % regressor 1
            allTrs(m,2) = 1;
            break;
        elseif timingData.preChoiceOff(n) <= tmp && tmp < timingData.infChoiceTime(n)
            % regressor 2
            allTrs(m,3) = 1;
            break;
        elseif timingData.drawBegTime(n) <= tmp && tmp < timingData.infoOn(n)
            % regressor 3
            allTrs(m,4) = 1;
            break;
        elseif timingData.freeDrawOn(n) <= tmp && tmp < timingData.freeDrawOff(n)
            % regressor 4
            allTrs(m,5) = 1;
            break;
        elseif timingData.betOn(n) <= tmp && tmp < timingData.fdbkOn(n)
            % regressor 5
            allTrs(m,6) = 1;
            break;
        elseif timingData.fdbkOn(n) <= tmp && tmp <= timingData.fdbkOff(n)
            % regressor 6
            allTrs(m,7) = 1;
            break;
        elseif timingData.tooSlowOnTime(n) <= tmp && tmp <= timingData.tooSlowOffTime(n)
            % regressor 7
            allTrs(m,8) = 1;
            break;
        else
        end
    end
end


fslregressor = struct;
infoData=getInfoValueFromBlock(trialData);
drawornot = trialData.draw-nanmean(trialData.draw);
leftorright = trialData.leftrightbet-nanmean(trialData.leftrightbet);
% choiceNames={'modTrialInfValue' 'modExpMutInfGain' 'modTrialStateValue'};
% infoNames={'modDrawActualMutInfoGain' 'modDrawActualValueGain' 'modDrawActualKLdivergence'};

choicePhaseVars=zscore([infoData.modTrialInfValue infoData.modExpMutInfGain infoData.modTrialStateValue]);
infoPhaseVars=[infoData.modDrawActualMutInfoGain infoData.modDrawActualValueGain infoData.modDrawActualKLdivergence];
choicePhaseVars=choicePhaseVars-repmat(nanmean(choicePhaseVars), length(choicePhaseVars), 1); 

% not sure if this code is getting used at all...will delete if unnecessary
% if ~unique(trialData.hiValSide)
%     beadDifference=trialData.startLeft-trialData.startRight;
% else
%     beadDifference=trialData.startRight-trialData.startLeft;
% end
    
% ll=length(infoData.modTrialInfValue);
% t=0; 

% allTrialDraws=nan(ll, 1);
infoOn=timingData.infoOn; 
infoTrials=isfinite(infoOn);
infoPhaseVars(infoTrials,:)=(infoPhaseVars(infoTrials,:)-repmat(nanmean(infoPhaseVars(infoTrials,:)), sum(infoTrials), 1) );

fslregressor.reg1 = struct;
fslregressor.reg1.modTrialInfValue = [timingData.preChoiceOn timingData.preChoiceOff-timingData.preChoiceOn choicePhaseVars(:,1)];
fslregressor.reg1.modExpMutInfGain = [timingData.preChoiceOn timingData.preChoiceOff-timingData.preChoiceOn choicePhaseVars(:,2)];
fslregressor.reg1.modTrialStateValue = [timingData.preChoiceOn timingData.preChoiceOff-timingData.preChoiceOn choicePhaseVars(:,3)];
fslregressor.reg1.drawBead = [timingData.preChoiceOn timingData.preChoiceOff-timingData.preChoiceOn drawornot];
fslregressor.reg1.box = [timingData.preChoiceOn timingData.preChoiceOff-timingData.preChoiceOn ones(length(choicePhaseVars(:,3)),1)];

fslregressor.reg2 = struct;
fslregressor.reg2.box = [timingData.preChoiceOff timingData.infChoiceTime-timingData.preChoiceOff ones(length(timingData.infChoiceTime),1)];

fslregressor.reg3 = struct;
fslregressor.reg3.modDrawActualMutInfoGain = [timingData.drawBegTime timingData.infoOn-timingData.drawBegTime infoPhaseVars(:,1)];
fslregressor.reg3.modDrawActualValueGain = [timingData.drawBegTime timingData.infoOn-timingData.drawBegTime infoPhaseVars(:,2)];
fslregressor.reg3.modDrawActualKLdivergence = [timingData.drawBegTime timingData.infoOn-timingData.drawBegTime infoPhaseVars(:,3)];
fslregressor.reg3.box = [timingData.drawBegTime timingData.infoOn-timingData.drawBegTime ones(length(infoPhaseVars(:,3)),1)];

fslregressor.reg3.modDrawActualMutInfoGain = fslregressor.reg3.modDrawActualMutInfoGain(~any(isnan(fslregressor.reg3.modDrawActualMutInfoGain),2),:);
fslregressor.reg3.modDrawActualValueGain = fslregressor.reg3.modDrawActualValueGain(~any(isnan(fslregressor.reg3.modDrawActualValueGain),2),:);
fslregressor.reg3.modDrawActualKLdivergence = fslregressor.reg3.modDrawActualKLdivergence(~any(isnan(fslregressor.reg3.modDrawActualKLdivergence),2),:);
fslregressor.reg3.box = fslregressor.reg3.box(~any(isnan(fslregressor.reg3.box),2),:);

fslregressor.reg4 = struct;
fslregressor.reg4.box = [timingData.freeDrawOn timingData.freeDrawOff-timingData.freeDrawOn ones(length(timingData.freeDrawOn),1)];
fslregressor.reg4.box = fslregressor.reg4.box(~any(isnan(fslregressor.reg4.box),2),:);

fslregressor.reg5 = struct;
fslregressor.reg5.box = [timingData.betOn timingData.fdbkOn-timingData.betOn ones(length(timingData.betOn),1)];
fslregressor.reg5.box = fslregressor.reg5.box(~any(isnan(fslregressor.reg5.box),2),:);
fslregressor.reg5.betleftright = [timingData.betOn timingData.fdbkOn-timingData.betOn leftorright];
fslregressor.reg5.betleftright = fslregressor.reg5.betleftright(~any(isnan(fslregressor.reg5.betleftright),2),:);

fslregressor.reg6 = struct;
fslregressor.reg6.box = [timingData.fdbkOn dataStruct.fdbkOff-timingData.fdbkOn ones(length(timingData.fdbkOn),1)];
fslregressor.reg6.box = fslregressor.reg6.box(~any(isnan(fslregressor.reg6.box),2),:);

fslregressor.reg7 = struct;
fslregressor.reg7.box = [timingData.tooSlowOnTime timingData.tooSlowOffTime-timingData.tooSlowOnTime ones(length(timingData.tooSlowOnTime),1)];
fslregressor.reg7.box = fslregressor.reg7.box(~any(isnan(fslregressor.reg7.box),2),:);
end