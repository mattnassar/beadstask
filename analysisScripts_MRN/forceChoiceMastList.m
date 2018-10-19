%% Code for Beads task behavioral analysis and EV file creation.

% 4-18-16

% store block info so we can ditch the last trial in each block. 


% 3-22-16

% added a section to look for univariate effects of house/face encoding using trial info coefficients. 
% first pass analysis looked terrible... realized that some of the
% coefficients are poorly estimated. you can see that by looking at the
% variance across coefficients estimated for different voxels for a single
% trial. High variance is indicative that there wasn't much to go on... I
% analyzed the data using two strategies... 1=remove high variance
% trials... 2= use a weighted regression with the weights computed
% according to proxy variance. First strategy revealed some blobs, but the
% motor asymmetry was not very clear... need to run the second strategy
% again because it didn't finish.


% 5-19-15

% changed directories for data... as apparently the data has been dumped in
% a different directory. 


% 2-16-15
% setup a new script for creating EVs that include a couple of real pretty
% obvious motor and visual controls (writeEVfile_fullModel.m)




%% get everything we need on the path.
clear classes

rootDir='/Users/mattnassar/matt_work_stuff/Matt/m_files/beadstask/analysisScripts_MRN'
cd(rootDir);
addpath(genpath(rootDir))

addpath(genpath(rootDir))

infoPath='/Users/mattnassar/Dropbox/BeadsTaskCode/EVs/infoTimes/'
%% load subject data
subFileDir='/Users/mattnassar/Dropbox/beadstaskcode/BeadsTask5data';
figDir='/Users/mattnassar/Dropbox/BeadsTaskCode/figures/'
trialCoeffDir= '~/Dropbox/BeadsTaskCode/trialCoeffs/';


subNames={'CC528', 'FL1136', 'FT3594', 'IA3593', 'IL3520',	'LL3555', ...
    'NS950', 'OF3592', 'TK3556', 'TQ3543', 'TQ3600', 'BA3550', 'KQ3548', 'MF3567',  'NN3554' };





% loop through each subject and load data, store in "allSubjData" structure
% and save that structure so we don't need to use this slow loading code
% again.


clear  allSimpData
for i = 1:length(subNames)
    fileName=fullfile(subFileDir, subNames{i});
    [allData]=beadTaskDataLoader_mrn(fileName);
    
    % put data structures in a bigger structure including all subjects.
    eval(sprintf('allSubjData.%s=allData;', subNames{i}));
    
    % beadTaskDataLoader changes directories, so lets go back manually.
    cd(subFileDir)


    
    
    
    dat1=unpackBeadsData(allData.realScannerGame1.block1.statusData);
    dat1.blk=ones(size(dat1.trialNum));
    dat1=straightStruct(dat1); 
    
    dat2=unpackBeadsData(allData.realScannerGame1.block2.statusData);
    dat2.blk=ones(size(dat2.trialNum)).*2;
    dat2=straightStruct(dat2);
 
    dat3=unpackBeadsData(allData.realScannerGame2.block1.statusData);
    dat3.blk=ones(size(dat3.trialNum)).*3;
    dat3=straightStruct(dat3);

    dat4=unpackBeadsData(allData.realScannerGame2.block2.statusData);
    dat4.blk=ones(size(dat4.trialNum)).*4;
    dat4=straightStruct(dat4);


    simpData=catBehav(dat2, dat1, true);
    simpData=catBehav(dat3, simpData, true);
    simpData=catBehav(dat4, simpData, true);       
    simpData.subj=i.*ones(size(simpData.trialNum));

    if ~exist('allSimpData')
        allSimpData=simpData;
    else
        allSimpData=catBehav(simpData, allSimpData, true);
    end
end


cd(figDir)


TR=1.5
TR_perBlock=395;
blkTime=TR.*TR_perBlock;

% Get data in favorite format?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%optionally, run commented lines and start from here:  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load allSubjData.mat
% subNames=fieldnames(allSubjData);


exSub='IL3520';
eval(sprintf('gameNames=fieldnames(allSubjData.%s);', exSub));
gameCat=1:length(gameNames);
allBlocks=fieldnames(allSubjData.IL3520.realScannerGame1);
allNames=fieldnames(allSubjData.IL3520.realScannerGame1.block1.statusData);



% Timing reminders:

% NOTE: dataStruct timings are relative to task start, but timingData
% timings are set relative to the start of the scanner.

% preChoiceOn    Subject can see the info choice. can't see mapping.
% preChoiceOff   Subject sees mapping and can make choice.

% infChoiceTime  Subject makes choice about whether he wants info.

% realInfoOn     Info is displayed in center of screen

% infoOn         counterintuitively, this actually encodes info ending

% freeDrawOn     begin free draw phase
% freeDrawOff    end free draw phase

% betOn          begin bet phase
% fdbkOn         end bet phase and begin feedback phase (no lag/jitter)

% fdbkOff        end feedback phase.


% So... data.infChoiceTime-data.preChoiceOff should be RT?



%% load all scanner data
scannedSubjs=subNames

%scannedSubjs={'HU2138', 'LX1803', 'NS969', 'OK547', 'QB2027', 'TA3514'};


blkNames={'realScannerGame1.block1', 'realScannerGame1.block2', 'realScannerGame2.block1', 'realScannerGame2.block2'};
evFolder='/Users/mattnassar/Dropbox/BeadsTaskCode/EVs/fullModel_EVs_2';
justHouseEvFolder='/Users/mattnassar/Dropbox/BeadsTaskCode/EVs/justHouse_EVs';
infoCoeffEvFolder='/Users/mattnassar/Dropbox/BeadsTaskCode/EVs/trialInfoEVs_EVs';


mkdir(justHouseEvFolder)
mkdir(evFolder)

% scan protocol info:
TR=1.5;
numTRs=400;


% loop through all subjects. Make EV files for each run, for each subject.

for i = 1:length(scannedSubjs)
    for j = 1:length(blkNames)
        
        % unpack data.
        eval(sprintf('data=unpackBeadsData(allSubjData.%s.%s.statusData);', ...
            scannedSubjs{i},  blkNames{j}));
        
        
        
        
        % get model based trial regressors.
        [numSubjs, uniqueVar, timingData, trialData, regData, allRegModVars]...
            =checkBlockEfficiency(data);
        data.infRT=data.infChoiceTime-data.preChoiceOff; % this is no longer relevant...
        close all
        
        % compute bead diff (high val beads - low value beads)
        if unique(trialData.hiValSide)==0
            beadDiff{i,j}=trialData.startLeft-trialData.startRight;
        else
            beadDiff{i,j}=trialData.startRight-trialData.startLeft;
        end
        
        
        % ok, i'm going to remove the "expected mutual information gain" as
        % i think that quantitiy is correlated with the actual kl
        % divergence
        
        allRegModVars.choice.vals= [data.curr_start_tokensRight+data.curr_start_tokensLeft, ...
            data.curr_start_tokensRight-data.curr_start_tokensLeft];
        allRegModVars.choice.vals=meanCentX(allRegModVars.choice.vals);
        allRegModVars.choice.names={'beadTotal', 'beadDifference'};
        
        allRegModVars.infOnRight=trialData.infoButtonSide;
        
        %         % make EV file
        %         evFileName=fullfile(evFolder, [scannedSubjs{i} '_EV_run' num2str(j) '_']);
        %         [regNames, timings]=writeEVfile(timingData, allRegModVars, evFileName);
        
        % make EV file
        evFileName=fullfile(evFolder, [scannedSubjs{i} '_EV_run' num2str(j) '_']);
        [regNames, timings]=writeEVfile_fullModel(timingData, allRegModVars, evFileName);
        
        %        keyboard
        
        mkdir(infoCoeffEvFolder);
        evFileName=fullfile(infoCoeffEvFolder, [scannedSubjs{i} '_EV_run' num2str(j) '_']);
        [regNames_trialCoef, timings, info]=writeEVfile_trialInfoCoeffs(timingData, allRegModVars, evFileName, trialData.trialNum);
        
        info.blk=repmat(j, size(info.on));
        info.on=info.on+ (info.blk-1).*blkTime;
        if j == 1
            subInfo=info;
        else
            %  keyboard
            subInfo=catBehav(info, subInfo, true);
        end
        
        
        %         evFileName=fullfile(justHouseEvFolder, [scannedSubjs{i} '_EV_run' num2str(j) '_']);
        %         [justHouse_regNames]=writeEVfile_justHouse(timingData, allRegModVars, evFileName);
        
        
        % add some modulators that are useful for our regression:
        % info on right = +1, info on left = -1.
        % took info     = +1, bet          = -1
        allRegModVars.infButtonPress.vals=  ((trialData.draw.*2)-1).*((trialData.infoButtonSide.*2)-1);
        allRegModVars.infButtonPress.names={'pushedRight'};
        
        
        allTimings{i,j}=timings;
        allModulators{i,j}=allRegModVars;
                
    end
    allInfo{i}=subInfo;
    
    % write a single file containing all Info events
    
    
    % Create a text file with all info times:
    infoTxt=subInfo.on
    
    
    fnInfo=fullfile(infoPath, sprintf('info_times_%s.txt', scannedSubjs{i}));
    dlmwrite(fnInfo, infoTxt, ' ');
    
    trialEVs=makeInfoTrialEVs(subInfo.on, TR, numTRs, subInfo.blk, scannedSubjs{i}, infoPath)

end

save beadsData.mat allTimings allModulators allInfo scannedSubjs



disp(regNames');
disp('done making EVs etc')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Whole brain GLM to identify house/face encoding.           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Goal: 1) identify regions that encode/house face information. First pass
% = univariate... just want to make sure we get FFA/PPA/ motor cortex

% does it make sense to do this in a univariate way? I didn't smooth all
% that much. I would still expect some effects (e.g. lateralized motor
% activity)... but i guess it is possible that everything lives in the
% voxel-by-voxel differences. 


% find last trial in each block:


isLastTrial=false(size(allSimpData.blk));
newBlock=false(size(allSimpData.trialNum));
for i = 1:length(allSimpData.blk)
    if i==1 || allSimpData.trialNum(i)<allSimpData.trialNum(i-1)
        newBlock(i)=true;
    else
        newBlock(i)=false;
    end
end
isLastTrial([find(newBlock(2:end)); length(newBlock)])=true;


%[isLastTrial, newBlock]



% issue is not just with the last trial. 2nd trial of 3rd block is crazy
% bad for subject number 1. Not sure what is going on there... it seems a
% bit odd that this particular trial would be so poorly estimated. 

% ditch poorly estimated trial coefficients?
varThresh=3;


a=dir([trialCoeffDir, '/*brainInfoTrialCoeffs.mat'])
fns={a.name};


firstInfo=nan(size(allSimpData.curr_trialInfList));
for i = 1:length(allSimpData.curr_trialInfList)
    if ~isempty(allSimpData.curr_trialInfList{i})
        firstInfo(i)=allSimpData.curr_trialInfList{i}(1);
    end
end

for i = 1:length(subNames)
    abNames{i}=subNames{i}(1:5);
end


infoCoeff=nan(91, 109,  91, length(fns)) ;
allSubTrialVar=nan(size(firstInfo));


for i = 1:length(fns)
     disp(sprintf('loading file: %s', fns{i}))       
    trialInfoCoeffs=load(fullfile(trialCoeffDir,fns{i}));
    % OK, have coefficients. Probably should mask with something?
    dims=size(trialInfoCoeffs.brainTrialCoeffs);
    selSubj=allSimpData.subj== find(strcmp(fns{i}(1:5),   abNames));

    
    % use whole brain response as a proxy for the variance on our estimate
    % of regressors (eg... if variance across brain is crazy, the the thing
    % is NOT well estimated:
    clear trialVar
    for k = 1:size(trialInfoCoeffs.brainTrialCoeffs, 4)
        trialSnap=trialInfoCoeffs.brainTrialCoeffs(:,:,:,k);
        trialVar(k)=nanstd(trialSnap(:));
    end
    
    dynThresh=min([prctile(trialVar, 96), varThresh]);
    sel= selSubj& isfinite(firstInfo);
%     hold on
%     plot(allSimpData.trialNum(selSubj& isfinite(firstInfo)))
%     plot(trialVar, 'r')

    trialVar(isLastTrial(selSubj& isfinite(firstInfo)))
    allSubTrialVar(selSubj&  isfinite(firstInfo))=trialVar;
    xes=meanCentX([ones(size(firstInfo(sel&allSubTrialVar<dynThresh))), firstInfo(sel&allSubTrialVar<dynThresh)]);
    tic
    for x= 1:dims(1)
        x
        for y=1:dims(2)
            for z=1:dims(3)
                dat=zscore(squeeze(trialInfoCoeffs.brainTrialCoeffs(x,y,z,trialVar<dynThresh)));
                if nanstd(dat)>0
                    %  B = regressW_mike(dat, trialVar',xes);
                    B = regress(dat, xes);
                    infoCoeff(x,y,z,i)=B(2);
                end
            end
        end
    end
    toc
end


% look at example subject
exCoeff=infoCoeff(:,:,:,14);
hf_sub_nii=make_nii(exCoeff);
view_nii(hf_sub_nii)
clear trialInfoCoeffs


tMap=nanmean(infoCoeff, 4)./(nanstd(infoCoeff,[],4)./sqrt(length(fns)));
varMap=nanstd(infoCoeff,[],4);
meanMap=nanmean(infoCoeff,4);

thresh=.01;
lb=icdf('t', thresh, 14);
ub=icdf('t', 1-thresh, 14);

    
% OK, initial results looked like crap. But that seems to be because there
% are some very poorly estimated trials. 

hf_nii=make_nii(squeeze(tMap), 2);

% save houseFaceRegression_4-19-16.mat


save_nii(hf_nii, 'houseFaceMap.nii')
%hf_nii.img(hf_nii.img>lb& hf_nii.img<ub)=0;

scale_t=hf_nii;
scale_t.img(scale_t.img>4)=4;
scale_t.img(scale_t.img<-4)=-4;

view_nii(scale_t)



houseFaceMap=load_nii('/Users/mattnassar/Dropbox/BeadsTaskCode/BeadsTask5data/houseFaceMap.nii')

%%%%%%%%%%%%%%%%%%%%% ROI Analysis %%%%%%%%%%%%%%%%%%%%%

% here i'm going to look for the level of "faciness" for each 
% information epoch. We'll do this by creating an 
% xmatrix that that has a separate regressor for every "info" phase


% notes from sunday 12-16-14:
% difference signals provide a reasonable level of separation between
% house and face info distributions. 
%
% however, the modulation of "discriminibility scores seems less
% compelling.

% there are several problems with the current analysis that could easily be
% corrected:

% 1) include infoChoice button time(and button pushed mod)
% 2) include subsequent button presses (ie extra info draws) in model
% 3) should pull main effects of inforational terms out of info series
% 4) outlier removal before classification? it seems like these outliers
% hurt



%% LOAD ROI data.


load allBeadsTseries_2-26-14.mat 
%load allBeadsTseries.mat
%load allBeadsTseries_old.mat


%% Create explanatory regressors.
% first pass we'll include:
% infoChoiceTime
% bet time
% too slow time
% bet time
% bet choice

% Now i'm going to add
% infoChoiceMade
% infoChoice --> pushed left





%% RUN REGRESSION FOR EACH SUBJECT.
%  use choice, bet, &motion parameters to soak up variability.
%  use individual regressors for each info event to get measure of
%  activation on each trial.

%  run regression model on data from each ROI.

% Plug in scan protocol info:
allScanTimes=TR.*.5+(0:TR:TR.*(numTRs-1));
numFixedRegs=13;  % should be 14 with infoChoice button push regressors?
int=.1;     % we'll model ongoing events with this level of precision. 
totRegThresh=10;
% first create block regressors will soak up overall level differences between
onBlock=ones(numTRs, 1);
offBlock=zeros(numTRs, 1);
blkRegs=[[onBlock; offBlock; offBlock; offBlock], ...
    [offBlock; onBlock; offBlock; offBlock], ...
    [offBlock; offBlock; onBlock; offBlock], ...
    [offBlock; offBlock; offBlock; onBlock]];

% these variables will be built over time, growing with each new subject
allBeadDiff   =[];
allDiscScores =[];
subNum        =[];
allZDiscScores=[];
allInfMods    =[];
allChoiceMods =[];
clear inSampAcc invP unbiasedAccMeasure 
for i = 1:length(scannedSubjs)-1
    sub_roiData=nans(numTRs.*4, 22);
    fixedRegs=nans(numTRs.*4, numFixedRegs);
    infRegs=nans(numTRs.*4, 0);
    
    % these variabls grow with each block and are reset for each new
    % subject
    infBeadDiff =[];
    subInfMods  =[]; %nan(0, 4);
    subInfRegs   =[];
    subChoiceMods=[];
    

    for j = 1:length(blkNames)
        
        eval(sprintf('run_roiData=allTseries.%s_r0%s;', scannedSubjs{i}, num2str(j)));
        eval(sprintf('run_confounds=allConfounds.%s_r0%s;', scannedSubjs{i}, num2str(j)));
        
        % get data necessary for regression. 
        choiceOn        =allTimings{i,j}(1).on;
        choiceDur       =allTimings{i,j}(1).duration;
        infoOn          =allTimings{i,j}(2).on;
        infoDur         =allTimings{i,j}(2).duration;
        betOn           =allTimings{i,j}(3).on;
        betDur          =allTimings{i,j}(3).duration;
        tooSlowOn       =allTimings{i,j}(4).on;
        tooSlowDuration =allTimings{i,j}(4).duration;
        choiceSide      =allModulators{i,j}.bet.vals;
        
%       % added by MRN 2-23-14
        infButtonPushOn   =allTimings{i,j}(5).on;
        infButtonPushDur  =allTimings{i,j}(5).duration;
        infPressSide      = allModulators{i,j}.infButtonPress.vals;
%       ok, this guy is not locked to trials, lets see how this goes...
        extraPushOn   =allTimings{i,j}(6).on;
        extraPushDur  =allTimings{i,j}(6).duration;




        clear totEffect
        % get choice timing regressor.
        lag=0:int:choiceDur;
        for z=1:length(lag)
                onTimes1=choiceOn+lag(z);
                effectSpread=nan(length(onTimes1), length(allScanTimes));
            for f = 1:length(onTimes1)
                effectSpread(f,:)=fast_fslgamma(allScanTimes-onTimes1(f),6,3);
            end
                totEffect(z,:)=nansum(effectSpread);
        end
        choiceReg=meanCentVar(nansum(totEffect));
        
%% problem... this timing variable seems to come on too late... OK, looks fixed...
        % choiceButtonPush... one time event
        clear totEffect
        onTimes1=infButtonPushOn;
        effectSpread=nan(length(onTimes1), length(allScanTimes));
        for f = 1:length(onTimes1)
            effectSpread(f,:)=fast_fslgamma(allScanTimes-onTimes1(f),6,3);
        end
        totEffect=nansum(effectSpread);
        modEffect=nansum(effectSpread.*repmat(infPressSide, 1, length(allScanTimes)));

        infChoicePushReg=meanCentVar(totEffect);
        infChoiceSideReg= meanCentVar(modEffect);
        
        % get bet timing regressor.
        clear totEffect
        lag=0:int:betDur;
        for z=1:length(lag)
            onTimes1=betOn+lag(z);
            effectSpread=nan(length(onTimes1), length(allScanTimes));
            for f = 1:length(onTimes1)
                effectSpread(f,:)=fast_fslgamma(allScanTimes-onTimes1(f),6,3);
            end
            
            totEffect(z,:)=nansum(effectSpread);
            modEffect(z,:)=nansum(effectSpread.*repmat(choiceSide, 1, length(allScanTimes)));
        end
        betReg=meanCentVar(sum(totEffect));
        choiceSideReg= meanCentVar(sum(modEffect));
        


%         % get extra info push timing regressor.
%         clear totEffect
%         onTimes1=extraPushOn;
%         effectSpread=nan(length(onTimes1), length(allScanTimes));
%         for f = 1:length(onTimes1)
%             effectSpread(f,:)=fast_fslgamma(allScanTimes-onTimes1(f),6,3);
%         end
%         totEffect=nansum(effectSpread);
%         extraInfChoicePushReg=zscore(totEffect);




        % get TOO SLOW timing regressor.
        clear totEffect
        lag=0:int:tooSlowDuration;
        for z=1:length(lag)
            onTimes1=tooSlowOn+lag(z);
            effectSpread=nan(length(onTimes1), length(allScanTimes));
            for f = 1:length(onTimes1)
                effectSpread(f,:)=fast_fslgamma(allScanTimes-onTimes1(f),6,3);
            end
            
            totEffect(z,:)=nansum(effectSpread);
        end
        tooSlowReg=meanCentVar(sum(totEffect));
        
        
%         %        Check variables:
%         %
%                 hold on
%                 plot(choiceReg, 'b')
%                 plot(betReg, 'r')
%                 plot(tooSlowReg, 'g')
%                 plot(choiceSideReg, 'm')
%         
        
        
        %% Now its time to get a fuckload of
        %  info regressors.  1 for each info presentation...
        
        clear totEffect
        lag=0:int:infoDur;
        sel=isfinite(infoOn);
        allInfoEffect=zeros(sum(sel), length(allScanTimes));
        
        for z=1:length(lag)
            onTimes1=infoOn(sel)+lag(z);
            effectSpread=nan(length(onTimes1), length(allScanTimes));
            for f = 1:length(onTimes1)
                effectSpread(f,:)=fast_fslgamma(allScanTimes-onTimes1(f),6,3);
            end
            allInfoEffect=allInfoEffect+effectSpread;
        end

%       %  CHECK TO see what regressors look like!
%       % EACH regressor should be a single event.
%       % the last such event usually occurs after scanning has stopped.
% 
%         hold on
%         plot(allInfoEffect(1,:))
%         plot(allInfoEffect(2,:), 'r')
%         plot(allInfoEffect(3,:), 'y')
%         plot(allInfoEffect(4,:), 'g')
%         plot(allInfoEffect(5,:), 'k')
%         plot(allInfoEffect(6,:), 'c')
%         plot(allInfoEffect(end-1,:), 'm')
%         plot(allInfoEffect(end,:), '--m')
        
        % we only want to keep trials where the response was captured
        % during scanning... usually the last one happens after the scanner
        % has stopped. 
        
        goodInf{i,j}=sum(allInfoEffect')>totRegThresh;
        allInfoEffect=meanCentX(allInfoEffect')';


        % put this block of fixed regressors and ROI timeseries' in the
        % correct block position across the full session.
        blkSel= 1+(j-1).*numTRs :(j).*numTRs;
        fixedRegs(blkSel, :)=[ones(numTRs, 1), choiceReg', betReg', tooSlowReg', choiceSideReg', infChoiceSideReg', infChoicePushReg', run_confounds];
        sub_roiData(blkSel,:)=run_roiData;



        % NOTE: this is the matrix of info timing regressors. 
        % each info phase gets its own regressor.
        % timing regressors are 0 for all times outside of the relevant
        % block. timing regressors grow in number as we add blocks, since
        % each new info event means a new timing regressor.
        newInfRegs=zeros(numTRs.*4, sum(goodInf{i,j}));
        newInfRegs(blkSel,:)=allInfoEffect(goodInf{i,j}, :)';
        subInfRegs=[subInfRegs newInfRegs];

        % we'll also want to pass along the bead differential associated
        % with each info phase. And house/face.

        % ok, now lets get all info phase modulator variables...
        % the selection array (sel) selects just info trial variables.
        infMods=allModulators{i,j}.info_raw.vals(sel,:);
        % and the goodInf variable selects just those where the scanner
        % data is collected
        gInfMods=infMods(goodInf{i,j}, :);
        % % and subInfMods puts it all in an array with length equal to the
        % number of goodInf trials for the entire session.
        subInfMods=[subInfMods; gInfMods];
        
        % ok, now lets get all choice phase modulators in the same way.
        choiceMods=allModulators{i,j}.choice_raw.vals(sel,:);
        gChoiceMods=choiceMods(goodInf{i,j}, :);
        subChoiceMods=[subChoiceMods; gChoiceMods];
        
        % follow same procedure with beadDiff
        infBeadDiffs=beadDiff{i,j}(sel);
        goodInfBeadDiffs=infBeadDiffs(goodInf{i,j});
        infBeadDiff=[infBeadDiff goodInfBeadDiffs'];
    end
    

    %% create 1 regression matrix for entire session
    fRegs=size([fixedRegs blkRegs], 2); % number of fixed regressors
    nInfRegs=size(infRegs, 2);          % number of trial info regressors
    % run regression for each subject.
    xmat=[fixedRegs blkRegs subInfRegs]; % create big x-matrix
    

       clear rpXMat
       for rp=1:size(subInfRegs,2)
            if rp==1
               fixInfReg=nanmean(subInfRegs(:,[rp+1:end]),2);
            elseif rp==size(subInfRegs,2)
               fixInfReg=nanmean(subInfRegs(:,[1:rp-1]),2);
            else
               fixInfReg=nanmean(subInfRegs(:,[1:rp-1, rp+1:end]),2);
            end
               trialInfReg=subInfRegs(:,rp);
               rpXMat(:,:,rp)=[fixedRegs blkRegs trialInfReg];
               
        end




    clear infBetas regBetas rpInfBeta
    for rr=1:size(sub_roiData, 2)
        [B,BINT,R,RINT]=regress(sub_roiData(:,rr), xmat); % run regression for each ROI
        infBetas(:,rr)=B(fRegs+1:end); % just store trial info coefficients for each ROI
        
    % yikes... seems like a lot of the coefficents are taking on a value of
    % 0.  better check into this:
    %    badList=find(B==0);
    % OK.  turns out that regress is noticing that columns of xmat are
    % linearly dependent and setting the values for dependent coefficients
    % to zero. This is bad. We'd like to be able to estimate EACH
    % coefficient. I think that the best way to do this is to take the russ
    % poldrack approach... ie run the regression once for each info trial.

        % this is a pretty slow way to do things. Ideally, the xMat's for
        % each regressor should be created out of the loop
        
        for rp=1:size(subInfRegs,2)
            [B,BINT,R,RINT]=regress(sub_roiData(:,rr), rpXMat(:,:,rp));
            if B(end)==0
            disp('problem with regression')
            keyboard
            end

            rpInfBeta(rp,rr)=B(end);
        end

        % argh. this should not be that difficult. Just trying to
        % regularize trial regressors. 
        [val, rank]=sort(rpInfBeta(rp,rr));     
        regBetas(rank,rr)=norminv((1:length(rank))./(max(rank)+1), 0, 1);          
    end

    house=subInfMods(:,4)==2;  % was the coin a house or face?
    hWeight= (house-.5).*-2; % houses are negative, faces are positive
    % this is for arrays that have data for all subjects
    ll=length(house);
    allBeadDiff(end+1:end+ll)=infBeadDiff;
    allInfMods(end+1:end+ll,:)=subInfMods;
    allChoiceMods(end+1:end+ll,:)=subChoiceMods;
    
    clear infBetaResid

    %xes=[ones(length(allInfMods),1), allInfMods(:,[1, 2, 3, 4]) allChoiceMods(:,[1 3])];
    % get rid of baseline variance for terms that we'll put in regression
    % model...
    xmat=[ones(size(house)), subInfMods(:,[2])];
    for rr=1:size(sub_roiData, 2)
        [~,~,infBetaResid(:,rr)]=regress(rpInfBeta(:,rr), xmat);
    end

   
    % now we are forcing normality... not sure if this is good.
    roiList={'ppa_ffa', 'ppa_ffa_wSel', 'motor', 'LR_ppa_ffa', 'LR_ppa_ffa_sel', 'allPreMotorEtc'};
    roiDat{1}=rpInfBeta(:,1:2);
    roiDat{2}=rpInfBeta(:,3:4);
    roiDat{3}=rpInfBeta(:,5:6);
    roiDat{4}=rpInfBeta(:,7:10);
    roiDat{5}=rpInfBeta(:,11:14);
    roiDat{6}=rpInfBeta(:,[5:6 15:22]);


   % keyboard
    clear discStat
    for rrr=1:length(roiList)
    [discStat(:,rrr), inSampAcc(rrr,i), invP(rrr,i) unbiasedAccMeasure(rrr,i) Betas]...
        =getDiscStat(house, roiDat{rrr}, 1000);
    end
    
    allDiscScores(end+1:end+ll,:)=(discStat);
    allZDiscScores(end+1:end+ll,:)=zscore(discStat)
    subNum(end+1:end+ll)=ones(1,ll).*i;
end
gSub=any((invP>.95));
disp('done creating discriminability index')



% % all data points...    
% boxplotM(allDiscScores, allBeadDiff',  .3, 1,  'r', 'k', [], [5 25 50 75 95]);

outlierThresh=3;
subNumThresh=2;
% plot means for each subject...
diffs=unique(allBeadDiff);
for i = 1:length(diffs)
binCount=histc(subNum(allBeadDiff==diffs(i)), 1:5);
dPoints(i)=sum(binCount>subNumThresh);
end
plot(diffs, dPoints);
close all
gDiffs=diffs(dPoints>4);



allMeans=nan(max(subNum), length(diffs));
allWeight=nan(max(subNum), length(diffs));


% create a regression matrix to explain discrinability scores
% NAMES IN : allRegModVars.info.names
%          : allRegModVars.choice.names

% combVar=zscore(allChoiceMods(:,1))+zscore(allChoiceMods(:,2))


%allNames={'intercept', allRegModVars.info.names{1:4}, allRegModVars.choice.names{[1,3]}} 
allRegModVars.info.names{5}='totInfo'
allNames={'intercept', allRegModVars.info.names{[5 2 4]}, allRegModVars.choice.names{[1,2, 3]}} 

xes=[ones(length(allInfMods),1), allInfMods(:,[1, 2,  4]) allChoiceMods(:,[1 2 3])];

% xes=[ones(length(allInfMods),1),  allChoiceMods(:,[3])];

surf2d(corr(xes))
set(gca, 'clim', [-1 1])
colorbar


clear B BINT allMeans allWeight resMeans
for i = 1:max(subNum)
    for k = 1:size(allDiscScores,2)
        subSel=subNum'==i & abs(allDiscScores(:,k))<outlierThresh;
        
        [B(i,:,k),BINT(i,:,:, k),R,RINT,STATS] = regress(allZDiscScores(subSel,k),meanCentX(xes(subSel,:)));
       

%             
%         [r, p]=  corr(allInfMods(subSel,3), allDiscScores(subSel,k));
%         plot(allInfMods(subSel,3), allDiscScores(subSel,k),  '.')
%         plot(allDiscScores(subSel,k))


        for j = 1:length(gDiffs)
            if j==length(gDiffs)
                tSel = allBeadDiff>=gDiffs(j);
            elseif j == 1
                tSel = allBeadDiff<=gDiffs(j);
            else
                tSel = allBeadDiff==gDiffs(j);
            end
            allMeans(i,j,k) = nanmean(allDiscScores(subSel&tSel', k));
            resMeans(i,j, k) = nanmean(R(tSel(subSel)));
            allWeight(i,j)= nansum(tSel'&subSel);
        end
        
    end
end
close all

realMeans=allMeans(:,:,4);
groupMean=nanmean(realMeans(gSub,:));
groupDev =nanstd(realMeans(gSub,:));
sel=sum(isfinite(realMeans))>1;


defaultPlotParameters
close all
hold on 
plot(gDiffs(sel), groupMean(sel), 'o', 'markerSize', 14, 'markerFaceColor', 'r', 'markerEdgeColor', 'k', 'lineWidth', 1)
%plot(gDiffs(sel), groupMeanMot(sel), 'o', 'markerSize', 14, 'markerFaceColor', 'g', 'markerEdgeColor', 'k', 'lineWidth', 1)


diffMat=repmat(gDiffs, max(subNum), 1);
plot(diffMat(:), realMeans(:),  'o', 'markerSize', 6, 'markerFaceColor', 'b', 'markerEdgeColor', 'k', 'lineWidth', 1)
xlim([-10 10])

% subplot(2, 1, 1)
 
hold on
for i = 1:size(realMeans, 1)
if(gSub(i))
plot(gDiffs, realMeans(i,:), 'color', cbColors(i,:))
end
end
plot(gDiffs(sel), groupMean(sel), 'o', 'markerSize', 14, 'markerFaceColor', 'r', 'markerEdgeColor', 'k', 'lineWidth', 1)


% previous workspace saved here:
%save beadsFMRI_workspace_2-23-14.mat
close all
inds=[5 6];
subYs=unbiasedAccMeasure(inds,:)' 
subSig=invP(inds, :)'
xVals=ones(size(subYs)).*     repmat(1:size(subYs, 2), size(subYs, 1), 1);
yMax= max(abs(subYs(:)));
yExtra= yMax./10
yLims=[-yMax-yExtra yMax+yExtra]
xLims=[min(xVals(:))-.5 max(xVals(:))+.5]
hold on
plot(xLims, [0 0], '--k')
plot(xVals, subYs, 'o', 'markerSize', 10, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
for i = 1:length(subSig(:))
if subSig(i)>.95
plot(xVals(i), subYs(i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(4,:), 'markerEdgeColor', 'k', 'lineWidth', 1); 
end
end
xlim(xLims)
ylim(yLims)
ylabel('Classifier accuracy score')
set(gca, 'box', 'off', 'XTickLabel', {})
saveas(gcf, 'classifierAccuracyFig.eps', 'epsc2')
close all



reg= 5
regName=allNames{reg}
inds=[5  6];
subYs=squeeze(B(:, reg, inds)) ;
[h p]=ttest(subYs)
subSig=invP(inds, :)';
xVals=ones(size(subYs)).*     repmat(1:size(subYs, 2), size(subYs, 1), 1);
yMax= max(abs(subYs(:)));
yExtra= yMax./10
yLims=[-yMax-yExtra yMax+yExtra];
xLims=[min(xVals(:))-.5 max(xVals(:))+.5];
hold on
plot(xLims, [0 0], '--k')
plot(xVals, subYs, 'o', 'markerSize', 10, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1);


% for i = 1:length(subSig(:))
% if subSig(i)>.95
% plot(xVals(i), subYs(i), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(4,:), 'markerEdgeColor', 'k', 'lineWidth', 1); 
% end
% end
title(regName);
xlim(xLims);
ylim(yLims);
ylabel('Coefficient')
set(gca, 'box', 'off', 'XTickLabel', {})
saveas(gcf, sprintf('accuracyMod_%s.eps', regName), 'epsc2');
close all








