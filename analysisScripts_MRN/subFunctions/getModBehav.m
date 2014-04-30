function [modBehav condInfoData, infoData, exValues]=getModBehav(dataStruct, invTemp)



ll=length(dataStruct.curr_maxUrnProb);
trialData.corrUrnSide=dataStruct.curr_urnType;              
trialData.maxUrnProb=dataStruct.curr_maxUrnProb;
trialData.startLeft=dataStruct.curr_start_tokensLeft;
trialData.startRight=dataStruct.curr_start_tokensRight;
trialData.infCost=dataStruct.curr_infCost;
trialData.hiValue=max([dataStruct.curr_rewCorrLeft dataStruct.curr_rewCorrRight], [], 2);
trialData.loValue=min([dataStruct.curr_rewCorrLeft dataStruct.curr_rewCorrRight], [], 2);
trialData.inValue=dataStruct.curr_penErrLeft;
trialData.curr_trialInfList=dataStruct.curr_trialInfList(1:ll);
trialData.draw=~cellfun(@isempty, trialData.curr_trialInfList);
trialData.hiValSide=dataStruct.curr_rewCorrRight>dataStruct.curr_rewCorrLeft;
trialData.betRight=dataStruct.curr_choice==2;

[infoData, modBehav, condInfoData, exValues]=getInfoValueFromBlock(trialData, invTemp);