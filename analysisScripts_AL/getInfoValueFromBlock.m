function infoData=getInfoValueFromBlock(trialData)
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
    infNet= drawValue-max(cat(3, pickLowVal, pickHighVal), [],3);
    infValue=infNet+infCost;
    sideBeads=[trialData.startLeft trialData.startRight];
    
    highValProbs=[1-q q];
    
    % now we loop only through trials for which the action value table is
    % relevant.
    for i = find(trialSel)'
        
        
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
                if isempty(trialData.curr_trialInfList{i})
                    keyboard
                end
                
                drawHigh=(trialData.hiValSide(i)+1)==trialData.curr_trialInfList{i}(1);
            else
                drawHigh=nan;
            end
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
            
            % if we draw toward the low value option
        else
            modDrawActualValueGain(i)=stateValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2)+trialData.inValue(i)-min(trialData.inValue)...
                -modTrialStateValue(i);
            modDrawActualMutInfoGain(i)=entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)-...
                entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2);
            
            pInit=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+1);
            pPost=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+2);
            modDrawActualKLdivergence(i)=computeKLdivInBits([pPost 1-pPost], [pInit 1-pInit]); 
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

