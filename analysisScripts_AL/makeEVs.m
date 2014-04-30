%% make change-point task EV's for FSL analysis.


evDir='/Users/mattnassar/matt_work_stuff/Matt/m_files/fMRI_analysis/localizerEVs/modEVs_3factor_10-11-13'
mkdir(evDir);
cd(evDir);

TR=2.5;
totTRs=226;
inFramesToDrop=5;
timeDropped=inFramesToDrop.*TR;


% need this to create nuisance variable text file.
numTRs=totTRs;
allScanTimes=0:TR:TR.*(numTRs-1);
noise=allDataStruct.currentStd;


%% CHANGES TO MAKE

%  1)  Set things up so we can easily drop a fixed number of trials from
%  the beginning of the session.
%  2)  Set up a way to add nuisance variables using the basic convolution.

%% MODEL:   outcome (nuisance), pCha, relUnc, outcomeValue 

% 2-4-13 model
%  modMat=[outcome modpCha modRelUnc vals LRres2];

% 2-13-13 model 
% GOAL: soak up all non-change-point related variance then look at the
% time series after onset of obvious change-points...

% 3-22-13 model
% GOAL: use model inspired by new behavioral analysis to dissociate LR
% effects. The behavioral model is set up in errVsPEbehavAnalysis.m

% 4-10-13 model
% GOAL: soak up non-change-point related variance and then look at the time
% series after onset of obvious changepoints.

% 5-14-13 model
% GOAL: soak up non-change-point related variance and then look at the time
% series after onset of obvious changepoints.  Use FSL 5.2.2.










%% find big / detectable change-points for timeseries analysis.
%  this model was for stick function regressors controlling time relative
%to detectable change-points
detCPs=modpCha>.5;
detCPs(newBlock)=false;
%modMat=[outcome vals detCPs];

bThresh=1;
cpMag=cat(1, nan, diff(allDataStruct.currentMean))./noise;
hist(abs(cpMag(CP)));
bCP=false(size(CP));
bCP((abs(cpMag)>bThresh))=true;
close all




% create a matrix of nuisance variables:

ll=length(prediction)

btLogic=false(ll,1);
btLogic(bt)=true;
pbtLogic=false(ll,1);
pbtLogic(bpt)=true;

jsLogic=false(ll,1);
jsLogic(joystickFail)=true;
pjsLogic=false(ll,1);
pjsLogic(joystickFail+1)=true;

nuisanceVars=double([btLogic pbtLogic jsLogic pjsLogic]);
%funcXes=[outcome (vals-.5) PE abs(CE)];


% identify trials where subject might have not fully executed joystick
% movement.

% 9-2-13 model
% GOAL: get activity pattern related to higher learning rates during
% outcome viewing.  Also model the residual updating (in direction of
% bucket)

%funcXes=[outcome abs(PE) predLR ySignedRes  PE abs(CE) nuisanceVars];

%% half of variance in LR explained by abs PE... 

clear rho

for i = 1:length(blockBegs)
    sel=zeros(length(funcXes), 1);
    sel( blockBegs(i):blockEnd(i))=true;
    gt=isfinite(predLR)&isfinite(ySignedRes)&sel;
    rho(:,:,i)=corr(funcXes(gt,:));
end


avRho=nanmean(rho, 3);

avRho.^2;

surf2d(avRho(1:end-4, 1:end-4))
colorbar
close all





% create a matrix of variables of interest

% MODEL BASED:
%errVsPEbehavAnalysis
%funcXes=[outcome modpCha  modRelUnc.*(1-modpCha) (vals-.5) PE abs(CE) nuisanceVars];

% model based with UNTRANSFORMED relative uncertainty
funcXes=[outcome modpCha  modRelUnc (vals-.5) PE CE.^2 nuisanceVars];



% model based JUST model predicted learning rates and binary residual term


%funcXes=[regModPredLR resDirection outcome PE abs(CE) nuisanceVars];










% ERROR BASED:  
% goal of this model is to distinguish signals relating change-point
% probability from those that simply reflect error magnitude, which is
% associated with a number of visual differences.

% the strategy is to use an outcome time regressor, then modulate this
% regressor according to a categorical variable (errTypeMatrix) that
% assigns each absolute error magnitude to a category.  The variable is
% exptrapolated such that errors falling between bin centers are assigned
% to both bins with weights determined by the relative distance to these
% bin centers.  errTypeMatrix is created in makeErrMagRegressors.m .

% MRN changed this... error regressor was colinear... now using really
% basic binning procedure (errTypeBin)... this definitely should not be
% colinear... 

% thresholds used for 6-26 EVs were:
% pTiles =
% (discarded) 2.3684    1) 5.5789    2) 9.2895   3) 13.6523   4) 19.1579   5) 28.2538   6) 47.5263  7) 277.0000
%

%funcXes=[outcome double(errTypeBin) (vals-.5) abs(CE) PE nuisanceVars];

imagesc(corr(funcXes))
colorbar
close all


% TIME BASED:
%funcXes=[outcome (vals-.5) PE abs(CE) nuisanceVars];


% Get rid of bad trials:

funcXes(allBT, :)=nan;
funcXes(allBT, 1)=outcome(allBT);
funcXes(~allDataStruct.isPredictionActive,:);
modMat=funcXes;

modulatorEVs=1;
timeSeriesEVs=0;
anyNuisance=0;
convNuisance=0;

% this loop will:
% 1: find subjects name... rename the subject if ID is too short
% 2: run through each block
% 3: select all data from that block
% 4: mean center each column of data
% 5: remove nan's for each column of data and replace with zero.



sub=struct


for i = 1:length(d)  % run through subjects
    
    % get subject ID for naming and selecting purposes
    subID=d(i).id;
    sel=strcmp(subID, allDataStruct.subjName); % select trials for a given subject 
    writeID=subID
    
    for j = 1:4 % run though task blocks
        blkNum=num2str(j); % block number (for naming purposes)
        sel1=strcmp(subID, allDataStruct.subjName)&allDataStruct.blkNum==j; % select trials for a given subject and block
        onTimes1=d(i).mrMain.outcomeTimes(sel1(sel));
        upTimes1=d(i).mrMain.updateTimes(sel1(sel));
        
        
        % get subject specific outcome duration;
        outcomeDuration=nanmean(upTimes1(2:end)-onTimes1(1:end-1));
        od(i)=outcomeDuration;
        updateDuration=nanmean(onTimes1-upTimes1);
        ud(i)=updateDuration;
      
        
        if j ==1
            sub.name{i}=subID;
            sub.firstStd(i)=unique(noise(sel1));
        end
        
        scaledEffect=nan(length(allScanTimes), length(onTimes1), size(nuisanceVars,2));
        if anyNuisance==1
            nVars=double(nuisanceVars(sel1,:));

            if any(any(nVars));
                convNVar=zeros(size(nVars));
                effectsSpread=nan(length(onTimes1), length(onTimes1));
                
                %% create a convolution matrix to show how an effect on one
                %% trial gets spread around
                
                clear effectSpread
                for f = 1:length(onTimes1)
                    effectTime=onTimes1(f);
                    effectSpread(:,f)=fast_fslgamma(allScanTimes-effectTime,6,3);
                    for n = 1:size(nVars,2)
                        scaledEffect(:,f,n) =   effectSpread(:,f).*nVars(i,n);
                    end
                end
            end
            
            nuisMat=squeeze(nansum(scaledEffect,2));
            nuisMat=meanCentX(nuisMat);
            
            if inFramesToDrop>0
                nuisMat=nuisMat(inFramesToDrop+1:end,:);
            end
            dlmwrite([writeID '_run' blkNum '_nuisance' ], nuisMat, ' ');
        end
        
    %% change onTimes1 to a different timeframe if we are going to drop
        %% frames
        
        
%         if any(nanstd(modMat(sel1,2:10))==0)
%             disp('houston: we have a problem')
%             keyboard
%         end
        
            
            
        
         
        for k=1:size(modMat, 2)+2
            if k==1
                featTxt=[onTimes1-timeDropped, ones(length(onTimes1),1), ones(length(onTimes1),1).*outcomeDuration];
            elseif k ==2
                featTxt=[upTimes1-timeDropped, ones(length(upTimes1),1), ones(length(onTimes1),1).*updateDuration];
            elseif modulatorEVs==1;
                selVals=modMat(sel1,k-2);   % grab data from the (k-1)th regressor
                selVals=selVals-nanmean(selVals); % mean center it
                selVals(~isfinite(selVals))=0; % replace nan's with 0;
                % and stick the values in FEAT's favorite format
                featTxt=[onTimes1-timeDropped, ones(length(onTimes1),1), ones(length(onTimes1),1).*selVals];
            end
                dlmwrite([writeID '_run' blkNum '_reg_' num2str(k)], featTxt, ' '); % write to a .txt file
        end
        
        
        if timeSeriesEVs==1;
            
            prevRegs=k;
            %% MAKE 3 column matrix!
            
            % find the times of big change-points for a given run
            timeRegs=onTimes1(bCP(sel1))-timeDropped;

            for kk=1:size(timeRegs,2);
                k=prevRegs+kk;
                timeRegK=timeRegs(:,kk);
                timeRegK=timeRegK(isfinite(timeRegK));
                featTxt=[timeRegK, ones(length(timeRegK),1).*TR, ones(length(timeRegK),1)];
                dlmwrite([writeID '_run' blkNum '_reg_' num2str(k)], featTxt, ' ');
            end 
        end
    end
end

