%getInfoValueFromModel


%% probability and likelihood of left and right
lLeft=log((trialData.maxUrnProb.^trialData.startLeft).*((1- trialData.maxUrnProb).^trialData.startRight));
lRight=log((trialData.maxUrnProb.^trialData.startRight).*((1- trialData.maxUrnProb).^trialData.startLeft));
pLeft=exp(lLeft-logsumexp([lLeft lRight], 2));
pRight=exp(lRight-logsumexp([lLeft lRight], 2));
sideProbs=[pLeft pRight];
hSide=unique(trialData.hiValSide+1);
lSide=mod(hSide, 2)+1;
expValChooseHigh=trialData.hiValue.*sideProbs(:,hSide);
expValChooseLow=trialData.loValue.*sideProbs(:,lSide);
tt=1;



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
        
        % here we will assume that the subject actually draws a bead
        
        if rand>(highValProbs(1+(trialData.corrUrnSide(i)==trialData.hiValSide(i)+1)))
            % bead favoring high value option drawn...
            modDrawActualValueGain(i)=stateValue(sideBeads(i,hSide)+2, sideBeads(i,lSide)+1)+trialData.inValue(i)-min(trialData.inValue)...
                -modTrialStateValue(i);
            modDrawActualMutInfoGain(i)=entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)-...
                entropy(sideBeads(i,hSide)+2, sideBeads(i,lSide)+1);
            
            
            pInit=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+1);
            pPost=pHighValUrn(sideBeads(i,hSide)+2,sideBeads(i,lSide)+1);
            
        else
            
            modDrawActualValueGain(i)=stateValue(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2)+trialData.inValue(i)-min(trialData.inValue)...
                -modTrialStateValue(i);
            modDrawActualMutInfoGain(i)=entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+1)-...
                entropy(sideBeads(i,hSide)+1, sideBeads(i,lSide)+2);
            
            pInit=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+1);
            pPost=pHighValUrn(sideBeads(i,hSide)+1,sideBeads(i,lSide)+2); 
        end
        modDrawActualKLdivergence(i)=computeKLdivInBits([pPost 1-pPost], [pInit 1-pInit]);
    end
end
    
    
 