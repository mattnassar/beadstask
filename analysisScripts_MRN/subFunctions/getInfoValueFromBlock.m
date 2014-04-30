function [infoData, modBehav, condInfoMetrics, exValues]=getInfoValueFromBlock(trialData, invTemp)
%% This function was originally written to get relevant informational quantities from
%  example task runs of the beads task. It now includes the option to spit
%  out model behavior for each trial. If you want model behavior, ask for 2
%  outputs and provide a softmax inverse temperature.
%  These additions were made by MRN 2/23/14

% modBehav is a structure with the following fields:
% choiceCategory   initial choice category: 1= pick low, 2= pick high, 3 =  draw.
% betSide          final bet side for initial draw choices.  1=low, 2=high.
% endDiff          findal difference in beads (high - low).

% If user asks for condInfoMetrics, compute them.
if nargout>2
    compConds=true;
else
    compConds=false;
end

if nargout>3
    getExVals=true;
else
    getExVals=false;
end




if nargin<2|isempty(invTemp)
    invTemp=10000;
end

if nargout<2
    runMod=0;
else
    runMod=1;
    modBehav.choiceCategory=nan(size(trialData.startRight));
    modBehav.betSide=nan(size(trialData.startRight));
    modBehav.endDiff=nan(size(trialData.startRight));
end



%% get info quantities for a given set of trials (in trialData)

%% trial data should have these fields.
% maxUrnProb
% startLeft
% startRight
% hiValSide (0 = left, 1 = right)
% infCost
% hiValue
% loValue
% inValue

%% optionally, if there is actual behavior attached, it can have these ones:
% curr_trialInfList
% draw

%getInfoValueFromModel
%keyboard

%% probability and likelihood of left and right
hSide=unique(trialData.hiValSide+1);
lSide=mod(hSide, 2)+1;

treeHeight=200;
q=unique(trialData.maxUrnProb);
infCost=unique(trialData.infCost);


%% compute an optimal action value table
%  THIS CODE modified by MRN in order to deal with blocks where there are
%  two sets of relevant conditions. In all cases the higher value should be
%  different, so we'll use that as a condition marker.

hCorrs=unique(trialData.hiValue);

modTrialStateValue=nan(size(trialData.startRight));
modTrialInfValue=nan(size(trialData.startRight));
modExpMutInfGain=nan(size(trialData.startRight));

modDrawActualMutInfoGain=nan(size(trialData.startRight));
modDrawActualValueGain=nan(size(trialData.startRight));
modDrawActualKLdivergence=nan(size(trialData.startRight));
infoType=nan(size(trialData.startRight));
stimInf=nan(size(trialData.startRight));


if compConds
    valGainGivHi=       nan(size(trialData.startRight));
    valGainGivLo=       nan(size(trialData.startRight));
    MutInfoGainGivHi=   nan(size(trialData.startRight));
    MutInfoGainGivLo=   nan(size(trialData.startRight));
    KLdivergenceGivHi=  nan(size(trialData.startRight));
    KLdivergenceGivLo=  nan(size(trialData.startRight));
    drawHiProb=         nan(size(trialData.startRight));
end


for j=1:length(hCorrs)
    hCorr=hCorrs(j);
    trialSel=trialData.hiValue==hCorr;
    
    % now get other variables JUST for this set of conditions..
    lCorr=unique(trialData.loValue(trialSel));
    hErr=unique(trialData.inValue(trialSel));
    lErr=unique(trialData.inValue(trialSel));
    
    
    [drawValue,  pickLowVal, pickHighVal, hThresh, lThresh, pHighValUrn, pHighValBallOnNextDraw] ...
        = btm_computeActionValueTable(treeHeight, q, infCost,  hCorr, lCorr,hErr, lErr);
    
    entropy=reshape(computeEntropy([pHighValUrn(:) 1-pHighValUrn(:)]), size(drawValue));
    [expEntropyAfterDraw]=computeFutureStateExpectation(entropy, pHighValBallOnNextDraw, 1-pHighValBallOnNextDraw);
    expLossInEntropy=entropy(1:end-1, 1:end-1)-expEntropyAfterDraw;
    stateValue= max(cat(3, pickLowVal, pickHighVal, drawValue), [],3);
    allVals= (cat(3, pickLowVal, pickHighVal, drawValue));
    infNet= drawValue-max(cat(3, pickLowVal, pickHighVal), [],3);
    infValue=infNet+infCost;
    sideBeads=[trialData.startLeft trialData.startRight];
    
    
    
    
    highValProbs=[1-q q];
    
    % now we loop only through trials for which the action value table is
    % relevant.
    for i = find(trialSel)'
        % if we want the model to do these trials, get its behavior
        if runMod
            % categories: 1= pick low, 2= pick high, 3 =  draw.
            hBeadsInd=sideBeads(i,hSide)+1;
            lBeadsInd=sideBeads(i,lSide)+1;
            relVals=squeeze(allVals(hBeadsInd, lBeadsInd, :));
            try
                pChoice=softMax(relVals, invTemp);
            catch
                keyboard
            end
            
            cumP=cumsum(pChoice);
            modBehav.choiceCategory(i)=find(rand<cumP, 1);
            % if model draws, lets see when it stops drawing
            % note: model draws from true urn probabilities... these
            % are not the same as the conditional bead probabilities.
            if modBehav.choiceCategory(i)==3
                probHigh=(highValProbs(1+(trialData.corrUrnSide(i)==trialData.hiValSide(i)+1)))
                modBehav.betSide(i)=nan;
                while ~isfinite(modBehav.betSide(i))
                    Draw=rand<probHigh;
                    if Draw
                        hBeadsInd=hBeadsInd+1;
                    else
                        lBeadsInd=lBeadsInd+1;
                    end
                    relVals=squeeze(allVals(hBeadsInd, lBeadsInd, :));
                    pChoice=softMax(relVals, invTemp);
                    cumP=cumsum(pChoice);
                    choice=find(rand<cumP, 1);
                    if choice<3  % if the choice is not to draw
                        modBehav.betSide(i)=choice;
                    end
                end
            end
            modBehav.endDiff(i)=hBeadsInd-lBeadsInd;
        end
        
        modTrialInfValue(i)=infValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1);
        modExpMutInfGain(i)=expLossInEntropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1);
        modTrialStateValue(i)=stateValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)+trialData.inValue(i)-min(trialData.inValue);
        
        if ~isfield(trialData, 'draw')
            % here we will assume that the subject actually draws a bead
            if rand>(highValProbs(1+(trialData.corrUrnSide(i)==trialData.hiValSide(i)+1)))
                % bead favoring high value option drawn...
                drawHigh=1;
            else
                drawHigh=0;
            end
        else
            if trialData.draw(i)
                
                infoType(i)=trialData.curr_trialInfList{i}(1);
                
                if isempty(trialData.curr_trialInfList{i})
                    keyboard
                end
                
                drawHigh=(trialData.hiValSide(i)+1)==trialData.curr_trialInfList{i}(1);
            else
                drawHigh=nan;
            end
        end
        
        % NOW compute the actual:
        %1) value gain. 2) MutInf gain and 3) KL divergence conditioned on draw:
        
        
        
        if compConds
            valGainGivHi(i)=stateValue(sideBeads(i,hSide)+2, sideBeads(i,lSide)+1)-stateValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1);
            valGainGivLo(i)=stateValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2)-stateValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1);
            MutInfoGainGivHi(i)=entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)-...
                entropy(sideBeads(i,hSide)+2, sideBeads(i,lSide)+1);
            MutInfoGainGivLo(i)=entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)-...
                entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2);
            pInit=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+1);
            pPost=pHighValUrn(sideBeads(i,hSide)+2,sideBeads(i,lSide)+1);
            KLdivergenceGivHi(i)=computeKLdivInBits([pPost 1-pPost], [pInit 1-pInit]);
            
            pInit=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+1);
            pPost=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+2);
            KLdivergenceGivLo(i)=computeKLdivInBits([pPost 1-pPost], [pInit 1-pInit]);
            
            drawHiProb(i)=pHighValBallOnNextDraw(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)
        end
        
        % conditionally, get example values for these trials.
        if getExVals
        exValues.pickLowVal(i)=pickLowVal(sideBeads(i,hSide), sideBeads(i,lSide));
        exValues.pickHighVal(i)=pickHighVal(sideBeads(i,hSide), sideBeads(i,lSide));
        exValues.drawValue(i)=drawValue(sideBeads(i,hSide), sideBeads(i,lSide));
        end

        
        
        
        
        % If there is no draw, set all quantities to nan
        if ~isfinite(drawHigh)
            modDrawActualValueGain(i)=nan;
            modDrawActualMutInfoGain(i)=nan;
            modDrawActualKLdivergence(i)=nan;
            % if we draw toward the high value option
        elseif drawHigh==1
            
            modDrawActualValueGain(i)=stateValue(sideBeads(i,hSide)+2, sideBeads(i,lSide)+1)+trialData.inValue(i)-min(trialData.inValue)...
                -modTrialStateValue(i);
            modDrawActualMutInfoGain(i)=entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)-...
                entropy(sideBeads(i,hSide)+2, sideBeads(i,lSide)+1);
            pInit=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+1);
            pPost=pHighValUrn(sideBeads(i,hSide)+2,sideBeads(i,lSide)+1);
            modDrawActualKLdivergence(i)=computeKLdivInBits([pPost 1-pPost], [pInit 1-pInit]);
            % added 2-18-14 by MRN
            stimInf(i)=log2(1./pHighValBallOnNextDraw(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1));
            % if we draw toward the low value option
        else
            
            modDrawActualValueGain(i)=stateValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2)+trialData.inValue(i)-min(trialData.inValue)...
                -modTrialStateValue(i);
            modDrawActualMutInfoGain(i)=entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)-...
                entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2);
            pInit=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+1);
            pPost=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+2);
            modDrawActualKLdivergence(i)=computeKLdivInBits([pPost 1-pPost], [pInit 1-pInit]);
            % added 2-18-14 by MRN
            stimInf(i)=log2(1./pHighValBallOnNextDraw(sideBeads(i,lSide)+1, sideBeads(i,hSide)+1));
            
        end
    end
end


% store expected info quantities.
infoData.modTrialInfValue=modTrialInfValue;
infoData.modExpMutInfGain=modExpMutInfGain;
infoData.modTrialStateValue=modTrialStateValue;

% store actual info quantities.
infoData.modDrawActualValueGain=modDrawActualValueGain;
infoData.modDrawActualMutInfoGain=modDrawActualMutInfoGain;
infoData.modDrawActualKLdivergence=modDrawActualKLdivergence;
infoData.infoType=infoType;
infoData.stimInf = stimInf;

if compConds
    condInfoMetrics.valGainGivHi=valGainGivHi;
    condInfoMetrics.valGainGivLo=valGainGivLo;
    condInfoMetrics.MutInfoGainGivHi=MutInfoGainGivHi;
    condInfoMetrics.MutInfoGainGivLo=MutInfoGainGivLo;
    condInfoMetrics.KLdivergenceGivHi=KLdivergenceGivHi;
    condInfoMetrics.KLdivergenceGivLo=KLdivergenceGivLo;
    condInfoMetrics.drawHiProb=drawHiProb;
end



