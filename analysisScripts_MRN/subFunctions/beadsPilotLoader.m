

load Minerva_realScannerGame1__1_blk201312191138beadsTaskData


%% TO FIX: trials definitely repeated in the end... 
%          Data saved after bad trials.
%          bad trials no longer replicated at end.

%          Delay issue: delay is close to 5... not sure if this is big
%          deal.

trialStartTime=[statusData.trialStartTime];
preChoiceOn=[statusData.preChoiceOn]
preChoiceOff=[statusData.preChoiceOff]

infChoiceOn=[statusData.infChoiceOn]
infChoiceMade=[statusData.infChoiceMade]
infoOff=[statusData.infoOff]

infoOn=[statusData.infoOn]
infoOff=[statusData.infoOff]

freeDrawOn=[statusData.freeDrawOn]
freeDrawOff=[statusData.freeDrawOff]

betOn=[statusData.betOn]

fdbkOn=[statusData.fdbkOn]
fdbkOff=[statusData.fdbkOff]


isGoodTrial=[statusData.isGoodTrial]
drawBegTime=[statusData.drawBegTime]



%% FIXED!


leftTok=[statusData.curr_start_tokensLeft]
rightTok=[statusData.curr_start_tokensRight]
leftRew=[statusData.curr_rewCorrLeft]

all=[leftTok' rightTok' leftRew']



