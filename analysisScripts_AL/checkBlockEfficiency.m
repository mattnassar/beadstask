%% TESTING THE QUALITY OF PARAMETERS DIRECTLY FROM FMRI BEADS TASK
clear classes
clc

load('NN3511_realScannerGame1__1_blk20141151759beadsTaskData.mat');
allN={statusData.name};
goodStatus=~cellfun(@isempty, allN);
statusData=statusData(goodStatus);
game11 = statusData;
load('NN3511_realScannerGame1__2_blk20141151814beadsTaskData.mat');
allN={statusData.name};
goodStatus=~cellfun(@isempty, allN);
statusData=statusData(goodStatus);
game12 = statusData;
load('NN3511_realScannerGame2__1_blk20141151827beadsTaskData.mat');
allN={statusData.name};
goodStatus=~cellfun(@isempty, allN);
statusData=statusData(goodStatus);
game21 = statusData;
load('NN3511_realScannerGame2__2_blk20141151838beadsTaskData.mat');
allN={statusData.name};
goodStatus=~cellfun(@isempty, allN);
statusData=statusData(goodStatus);
game22 = statusData;
% 
allNames=fieldnames(statusData);

%%
dataStruct.curr_trialInfList={statusData.curr_trialInfList};

for j = 1:length(allNames)
    % get data from a field
    if (isnumeric(eval(sprintf('statusData(1).%s', char(allNames(j)))))||...
            islogical(eval(sprintf('statusData(1).%s', char(allNames(j))))))...
            && length((eval(sprintf('statusData(1).%s', char(allNames(j))))))==1;
        eval(sprintf('dataStruct.%s=cat(1, statusData.%s);', char(allNames(j)), char(allNames(j))));
    elseif (isnumeric(eval(sprintf('statusData(1).%s', char(allNames(j)))))||...
            islogical(eval(sprintf('statusData(1).%s', char(allNames(j))))))...
            && length((eval(sprintf('statusData(1).%s', char(allNames(j))))))>1;
        eval(sprintf('dataStruct.%s={statusData.%s};', char(allNames(j)), char(allNames(j))));
    end
end

dataStruct=straightStruct(dataStruct);


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

inT=dataStruct.startTime;
timingData=struct;
timingData.infoOn=dataStruct.infoOn-inT;
timingData.choiceOn=dataStruct.preChoiceOn-inT;
timingData.feedbackOn=dataStruct.fdbkOn-inT;
timingData.infoOff=dataStruct.infoOff-inT;
timingData.choiceOff=dataStruct.preChoiceOff-inT;
timingData.feedbackOff=dataStruct.fdbkOff-inT;
timingData.trialStart=dataStruct.trialStartTime-inT;

%getInfoValueFromModel
[numSubjs, uniqueVar] = getRunEfficiency(trialData, timingData)


