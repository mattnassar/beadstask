%% load all scanner data
clear classes
cd /Users/mattnassar/Dropbox/BeadsTaskCode/beadsTask


codeDir='~/matt_work_stuff/Matt/m_files/beadstask/analysisScripts_MRN'
addpath(genpath(codeDir))



load allSubjData.mat
subNames=fieldnames(allSubjData);


behavSubjs={'BJ2376',	'BX2031',	'DT2005',	'NC3509',	'SQ3510','DZ572',	'HE3515',	'KQ3513',	'LT3516',	'LZ3512',	'NN3511', ...
    'OK547',	'SI3517',   'EO3518', 	'HU2138',	'LL2026',	'LX1803', ...
    'MC1906',	'NS969',	'QB2027',	'TA3514'};

noDrawBlkNames={'nodraw.block1', 'nodraw.block2'};
blkNames={'game1.block1', 'game1.block2', 'game2.block1', 'game2.block2'}; %
blk={'game1', 'game1', 'game2', 'game2'};

nonUndergrads={'LX1803', 'QB2027', 'HU2138', 'LL2026', 'MC1906', 'LT3516', 'BX2031', 'DT2005'};

%blkNames={'realScannerGame1.block1', 'realScannerGame1.block2', 'realScannerGame2.block1', 'realScannerGame2.block2'};
% blk={'realScannerGame1', 'realScannerGame1', 'realScannerGame2', 'realScannerGame2'}

makeFigs=false;
invT=1;

biasSum=nan(size(behavSubjs));
biasMax=nan(size(behavSubjs));
stopSum=nan(size(behavSubjs));

hiValStopSum  =nan(size(behavSubjs));
loValStopSum  =nan(size(behavSubjs));
hiValindPoint =nan(size(behavSubjs));
loValindPoint =nan(size(behavSubjs));
hiValDrawFrac =nan(size(behavSubjs));
loValDrawFrac =nan(size(behavSubjs));


totDrawFrac   =nan(size(behavSubjs));
totModDrawFrac=nan(size(behavSubjs));
totThreshDist =nan(size(behavSubjs));
totModThreshDist=nan(size(behavSubjs));



noDrawBias=nan(size(behavSubjs));
indPoint=nan(size(behavSubjs));
infoFitBias=nan(size(behavSubjs));
stopFitSum=nan(size(behavSubjs));
nonUndergradVar=nan(size(behavSubjs));
infoFitStd=nan(size(behavSubjs));
avThresh=nan(size(behavSubjs));
modInfoFitBias=nan(size(behavSubjs));
modStopSum=nan(size(behavSubjs));
defaultPlotParameters

for i = 1:length(behavSubjs)
    nonUndergradVar(i)=~isempty(strmatch(behavSubjs{i}, nonUndergrads));
    clear data allData noDrawData allNoDrawData
    for j = 1:length(noDrawBlkNames)
        % unpack data.
        eval(sprintf('noDrawData=unpackBeadsData(allSubjData.%s.%s.statusData);', ...
            behavSubjs{i},  noDrawBlkNames{j}));
        if ~isempty(noDrawData)
            noDrawData.hiValSide=(noDrawData.curr_rewCorrRight>noDrawData.curr_rewCorrLeft)+1;
            ll=length(noDrawData.trialNum);
            noDrawData.curr_trialInfList=noDrawData.curr_trialInfList(1:ll);
            noDrawData.anyDraw= ~cellfun(@isempty, noDrawData.curr_trialInfList);
            % compute bead diff (high val beads - low value beads)
            if unique(noDrawData.hiValSide)==1
                noDrawData.begBeadDiff=noDrawData.curr_start_tokensLeft-noDrawData.curr_start_tokensRight;
                noDrawData.endBeadDiff=noDrawData.curr_tokensLeft-noDrawData.curr_tokensRight;
            else
                noDrawData.begBeadDiff=noDrawData.curr_start_tokensRight-noDrawData.curr_start_tokensLeft;
                noDrawData.endBeadDiff=noDrawData.curr_tokensRight-noDrawData.curr_tokensLeft;
            end
            
            if j ==1
                allNoDrawData=straightStruct(noDrawData);
            else
                allNoDrawData=catBehav(noDrawData, allNoDrawData, 1);
            end
            
        end
    end
    
    noDrawBias(i)=sum((allNoDrawData.curr_choice==allNoDrawData.hiValSide))./length(allNoDrawData.hiValSide);
    %% it would be good to get a better measure of the indifference point.
    xes=allNoDrawData.begBeadDiff;
    yes=allNoDrawData.curr_choice==allNoDrawData.hiValSide;
    % Now we have one. Logistic regression to large reward choice data,
    % compute value of x where logistic function = 0 (ie probability = .5)
    B=glmfit(xes, yes, 'binomial');  % fit logistic function to data
    diffs=unique(allNoDrawData.begBeadDiff);
    modPred=B(1)+B(2).*diffs; % get logistic predictions from model fit
    indPoint(i)=-B(1)./B(2);  % calculate indifference point based on fit
    
    
    
    for j = 1:length(blkNames)
        % check whether data exists:
        if  eval(sprintf('isfield(allSubjData.%s, ''%s'');', ...
                behavSubjs{i},  blk{j}));
            % unpack data.
            eval(sprintf('data=unpackBeadsData(allSubjData.%s.%s.statusData);', ...
                behavSubjs{i},  blkNames{j}));
            
            if ~isempty(data)
                
                data.hiValSide=(data.curr_rewCorrRight>data.curr_rewCorrLeft)+1;
                
                ll=length(data.trialNum);
                data.curr_trialInfList=data.curr_trialInfList(1:ll);
                data.anyDraw= ~cellfun(@isempty, data.curr_trialInfList);
                
                % compute bead diff (high val beads - low value beads)
                if unique(data.hiValSide)==1
                    data.begBeadDiff=data.curr_start_tokensLeft-data.curr_start_tokensRight;
                    data.endBeadDiff=data.curr_tokensLeft-data.curr_tokensRight;
                else
                    data.begBeadDiff=data.curr_start_tokensRight-data.curr_start_tokensLeft;
                    data.endBeadDiff=data.curr_tokensRight-data.curr_tokensLeft;
                end
                
                [modData, condInfoTerms, infoData, exValues]=getModBehav(data, invT)
                
                % optionally, look at information metrics:
                if makeFigs
                    
                    
                    hold on
                    [orderedDiffs, ind]=sort(data.begBeadDiff);
                    wt=condInfoTerms.drawHiProb;
                    wHi=condInfoTerms.KLdivergenceGivHi;
                    wLo=condInfoTerms.KLdivergenceGivLo;
                    wDiff=wHi.*wt+wLo.*(1-wt);
                    ran=[min([wHi; wLo]) max([wHi; wLo])];
                    plot(orderedDiffs, wHi(ind), 'b')
                    plot(orderedDiffs, wLo(ind), 'g')
                    %            plot(orderedDiffs, wDiff(ind).*10, 'r')
                    plot([0 0], ran, '--k')
                    plot([-10 10], [0 0], '--k')
                    aa=legend('high val', 'low val');
                    set(aa, 'box', 'off', 'location', 'northeast')
                    xlabel('bead difference')
                    ylabel('KL Divergence')
                    set(gca, 'box', 'off')
                    saveas(gcf, 'KLDivergenceFig.eps', 'epsc2')
                    close all
                
                    hold on
                    [orderedDiffs, ind]=sort(data.begBeadDiff);
                     wDiff=wHi.*wt+wLo.*(1-wt);
                    ran=[0 1];
                    plot([0 0], ran, '--k')
                    a=plot(orderedDiffs, infoData.modTrialInfValue(ind)./max(infoData.modTrialInfValue), 'color', cbColors(4,:))
                    b=plot(orderedDiffs, infoData.modExpMutInfGain(ind)./max(infoData.modExpMutInfGain), 'color', cbColors(3,:))
                    %            plot(orderedDiffs, wDiff(ind).*10, 'r')
                    plot([-10 10], [0 0], '-k')
                    xlim([-10 10])
                    ylim([0 1])
                    aa=legend([a b], 'Info value', 'Info quantity');
                    set(aa, 'box', 'off', 'location', 'northeast')
                    xlabel('bead difference')
                    ylabel('Info value/quantity')
                    set(gca, 'box', 'off')
                    saveas(gcf, 'infoDissociationFig.eps', 'epsc2')
                    close all
      
                    hold on
                    sel=max([data.curr_rewCorrLeft data.curr_rewCorrRight]')==70
                    [orderedDiffs, ind]=sort(data.begBeadDiff(sel));
                           ran=[0 70];
                    plot([0 0], ran, '--k', 'lineWidth', 5)
                    
                    selExValues=selBehav(exValues, sel)
                    
                    c=plot(orderedDiffs, selExValues.drawValue(ind), 'color', 'g', 'lineWidth', 4)     
                    a=plot(orderedDiffs, selExValues.pickLowVal(ind), 'color', 'b', 'lineWidth', 4)
                    b=plot(orderedDiffs, selExValues.pickHighVal(ind), 'color', 'r', 'lineWidth', 4)

                   
                    xlim([-10 10])
                    ylim([ran])
                    aa=legend([a b c], 'Bet High', 'Bet Low', 'Draw');
                    set(aa, 'box', 'off', 'location', 'northwest')
                    xlabel('Bead difference')
                    ylabel('Value')
                    set(gca, 'box', 'off')
                    saveas(gcf, 'drawAndBetVals.eps', 'epsc2')
                    
            
                    
                    
                    
                    
                    
                    
               
                
                
                
                
                end
                
                
                
                if j ==1
                    allData=straightStruct(data);
                    allModData=straightStruct(modData);
                else
                    allData=catBehav(data, allData, 1);
                    allModData=catBehav(modData, allModData, 1);
                end
                
            end
        end
    end
    
    
    if exist('data', 'var')&&sum(allData.anyDraw)>0
        %% make classic plot.
        
        hiVal=max([allData.curr_rewCorrRight allData.curr_rewCorrRight]')';
        isHi =hiVal==max(hiVal);
        
        diffs=unique(allData.begBeadDiff);
        for k = 1:length(diffs);
            sel=allData.begBeadDiff==diffs(k);
            meanDrawFrac(k)=nanmean(allData.anyDraw(sel));
            meanBetHiFrac(k)=nanmean(~allData.anyDraw(sel)&(allData.hiValSide(sel)==allData.curr_choice(sel)));
            meanBetLoFrac(k)=nanmean(~allData.anyDraw(sel)&(allData.hiValSide(sel)~=allData.curr_choice(sel)));
            
            meanModDrawFrac(k)=nanmean(allModData.choiceCategory(sel)==3);
            meanModBetHiFrac(k)=nanmean(allModData.choiceCategory(sel)==2);
            meanModBetLoFrac(k)=nanmean(allModData.choiceCategory(sel)==1);
            
            
        end
        
        
        
        
        
        
        biasSum(i)=nanmean(allData.begBeadDiff(allData.anyDraw));
        biasMax(i)=nanmean(diffs(meanDrawFrac==max(meanDrawFrac)));
        totDrawFrac(i)=nanmean(allData.anyDraw);
        totModDrawFrac(i)=nanmean(allModData.choiceCategory==3);
        hiValDrawFrac(i)=nanmean(allData.anyDraw(isHi));
        loValDrawFrac(i)=nanmean(allData.anyDraw(~isHi));
        
        [estimates, sse] = fitGaussPDF_wGain (allData.begBeadDiff, allData.anyDraw, []);
        infoFitBias(i)=estimates(1);
        infoFitStd(i)=estimates(2);
        fitProb=estimates(4)+(estimates(3).*(1/sqrt(2*pi*estimates(2)) * exp(-(diffs-estimates(1)).^2/(2*estimates(2)))));
        
        [modEstimates, sse] = fitGaussPDF_wGain (allData.begBeadDiff, allModData.choiceCategory==3, []);
        modInfoFitBias(i)=modEstimates(1);
        modFitProb=modEstimates(4)+(modEstimates(3).*(1/sqrt(2*pi*modEstimates(2)) * exp(-(diffs-modEstimates(1)).^2/(2*modEstimates(2)))));
        
        
        % OK, now that we have a gaussian info distribution, lets use that
        % to identify "thresholds" in a way that is relative insensitive to
        % noise. 
  
        modTEsts(i,:)=getThreshFromGaussian(modEstimates)
        subTEsts(i,:)=getThreshFromGaussian(estimates)

  
  
        
        if makeFigs
            hold on
            plot([indPoint(i) indPoint(i)], [0 1], 'c')
            plot(diffs, meanDrawFrac, 'g')
            plot(diffs, meanBetHiFrac, 'r')
            plot(diffs, meanBetLoFrac, 'b')
            plot([0 0], [0 1], '--k')
            xlabel('Bead difference')
            ylabel('Choice frequency')
            plot(diffs, fitProb, '--m')
            
            set(gca, 'box', 'off')
            saveas(gcf, sprintf('indSubDrawPlot_%s.eps', behavSubjs{i}), 'epsc2')
            close all
        end
        
        %
        % plot model behavior?
        hold on
        plot([indPoint(i) indPoint(i)], [0 1], 'c')
        plot(diffs, meanModDrawFrac, 'g')
        plot(diffs, meanModBetHiFrac, 'r')
        plot(diffs, meanModBetLoFrac, 'b')
        plot([0 0], [0 1], '--k')
        plot(diffs, modFitProb, '--m')
        xlabel('Bead difference')
        ylabel('Choice frequency')
        set(gca, 'box', 'off')
        close all
        
        
        
        
        
        
        %% map thresholds based on end of betting:
        allPossDiffs=-30:30;
        sel=allData.anyDraw&(allData.hiValSide==allData.curr_choice);
        if sum(sel)>1
            hCount=histc(allData.endBeadDiff(sel), allPossDiffs);
            [estimates, sse] = fitGaussPDF_wGain (allPossDiffs', hCount, []);
            hDistFit=estimates(1);
        else
            hDistFit=nan;
        end
        
        
        sel2=allData.anyDraw&(allData.hiValSide~=allData.curr_choice);
        if sum(sel2)>1
            lCount=histc(allData.endBeadDiff(sel2), allPossDiffs);
            [estimates, sse] = fitGaussPDF_wGain (allPossDiffs', lCount, []);
            lDistFit=estimates(1);
        else
            lDistFit=nan;
        end
        
        stopFitSum(i)=hDistFit+lDistFit;
        
        
        
        highStop = nanmean(allData.endBeadDiff(allData.anyDraw&(allData.hiValSide==allData.curr_choice)));
        lowStop  = nanmean(allData.endBeadDiff(allData.anyDraw&(allData.hiValSide~=allData.curr_choice)));
        % compute same thing for high and low value trials separately
        hiHighStop = nanmean(allData.endBeadDiff(isHi&allData.anyDraw&(allData.hiValSide==allData.curr_choice)));
        loHighStop = nanmean(allData.endBeadDiff(~isHi&allData.anyDraw&(allData.hiValSide==allData.curr_choice)));
        hiLowStop  = nanmean(allData.endBeadDiff(isHi&allData.anyDraw&(allData.hiValSide~=allData.curr_choice)));
        loLowStop  = nanmean(allData.endBeadDiff(~isHi&allData.anyDraw&(allData.hiValSide~=allData.curr_choice)));
        
        modHighStop= nanmean(allModData.endDiff(allModData.choiceCategory==3&(allModData.betSide==2)));
        modLowStop= nanmean(allModData.endDiff(allModData.choiceCategory==3&(allModData.betSide==1)));
        




        totThreshDist(i)=highStop-lowStop;
        stopSum(i)=highStop+lowStop;  % should be zero if unbiased, negative if biased toward drawing in low val space
        avThresh(i)=highStop-lowStop;


        hiValStopSum(i)    =hiHighStop+hiLowStop;
        loValStopSum(i)    =loHighStop+loLowStop;
        modStopSum(i)      =modHighStop+modLowStop;
        totModThreshDist(i)=modHighStop-modLowStop;        


numBins=15; %max(allData.endBeadDiff)-min(allData.endBeadDiff);
        
        if makeFigs
            if (sum(sel2)>1)&&(sum(sel)>1)
                mCount=max([hCount; lCount]);
                
                hold on
                b=bar(allPossDiffs,lCount,'histc')
                set(b, 'faceColor', 'b')
                mExt=1+max(abs([min(allData.endBeadDiff) max(allData.endBeadDiff)]))
                xlim([-mExt mExt])
                plot([0 0], [0 mCount], '--k')
                xlabel('Bead difference')
                ylabel('frequency')
                set(gca, 'box', 'off')
                saveas(gcf, sprintf('indSubEndBeadHist_justLowVal_%s.eps', behavSubjs{i}), 'epsc2')
                
                a=bar(allPossDiffs,hCount,'histc')
                set(a, 'faceColor', 'r')
                
                saveas(gcf, sprintf('indSubEndBeadHist_%s.eps', behavSubjs{i}), 'epsc2')
                close all
            end
        end
    end
end


disp('done with basic anlalysis')

%%%%%%%%%%%%%%%%%%%%% GROUP ANALYSIS %%%%%%%%%%%%%%%%%%%%%

% quickie computation of optimal choice task behavior.
likeRatio=log((.6./.4).^diffs);
pHigh=1./(1+exp(-likeRatio))
evHigh=pHigh.*70;
evLow =(1-pHigh).*10;
optIndiff=nanmean([min(diffs(evHigh>evLow)) max(diffs(evHigh<evLow))])



% ideal:
% modStopSum = -7
% modInfoFitBias = -3.0187

optStopSum=-7;      % models stopping bias
optbiasSum=-3.0187; % this is actually the models fit bias

mixBias=nanmean([infoFitBias./optbiasSum; stopSum./optStopSum]);
modMixBias=nanmean([modInfoFitBias./optbiasSum; modStopSum./optStopSum]);



optMixBias=nanmean([optStopSum optbiasSum]);
sel=isfinite(stopSum)
hold on
plot([-4 4], [-4 4], '--k')


plot(stopFitSum(sel), stopSum(sel), '.b')
xlabel('fitStopSum')
ylabel('stopSum')
close all

infoFitBias(sel)




tt=5
hold on
plot([-tt tt], [0 0], '--k')
plot( [0 0], [-tt.*2 tt.*2], '--k')
plot(infoFitBias(sel), stopSum(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
%plot(modInfoFitBias(sel), modStopSum(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(6,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
%plot(modInfoFitBias(sel), modStopSum(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(8,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
ylabel('Stopping point bias')
xlabel('Initial draw bias')
set(gca, 'box', 'off')
saveas(gcf, 'groupBiasSummary.eps', 'epsc2')
plot(-3.0187, -7, 'o', 'markerSize', 14, 'markerFaceColor', 'r', 'markerEdgeColor', 'k', 'lineWidth', 1)
saveas(gcf, 'groupBiasSummary_wOpt.eps', 'epsc2')

plot(modInfoFitBias(sel), modStopSum(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(8,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
xlim([-5 5]);
ylim([-10 10]);



saveas(gcf, 'groupBiasSummary_wMod.eps', 'epsc2')
close all


close all
hold on
plot([-10 1], [0 0], '--k')
plot([0 0], [-1 1.5],  '--k')
% 
% plot([-2 -2], [-5 5], 'c')
nonUndergradVar=logical(nonUndergradVar);

[rho p]=corr(indPoint(sel)', mixBias(sel)')


plot(indPoint(sel), mixBias(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
%plot(indPoint(nonUndergradVar), mixBias(nonUndergradVar), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
ylabel('Information bias')
xlabel('Indifference point')
set(gca, 'box', 'off')
saveas(gcf, 'groupBiasVsIndiffPoint.eps', 'epsc2')
plot(-5, 1,  'o', 'markerSize', 14, 'markerFaceColor', 'r', 'markerEdgeColor', 'k', 'lineWidth', 1)
saveas(gcf, 'groupBiasVsIndiffPoint_wOpt.eps', 'epsc2')
close all






% 
% totDrawFrac   =nan(size(behavSubjs));
% totModDrawFrac=nan(size(behavSubjs));
% totThreshDist =nan(size(behavSubjs));
% totModThreshDist=nan(size(behavSubjs));


startThreshBias=sum(subTEsts, 2);
startThreshDist=(subTEsts(:,1)-subTEsts(:,2));

plot(startThreshBias, log(startThreshDist), '.')



hold on
plot([0 1], [0 0], '--k')
plot(totDrawFrac(sel), mixBias(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
set(gca, 'box', 'off')
ylabel('Information bias')
xlabel('Draw fraction')
saveas(gcf, 'drawsVsBias.eps', 'epsc2')
plot(totModDrawFrac(sel), modMixBias(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(8,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
saveas(gcf, 'drawsVsBias_wMod.eps', 'epsc2')
close all



hold on
plot(totDrawFrac(sel), totThreshDist(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1)
plot(totModDrawFrac(sel), totModThreshDist(sel), 'o', 'markerSize', 10, 'markerFaceColor', cbColors(8,:), 'markerEdgeColor', 'k', 'lineWidth', 1)




%% previous analysis workspace saved here:
%  save beadsBehavAnalysisWorkspace_2-13-14.mat


%% LEARNING EFFECTS ANALYSIS
%  LEARNING in noDraw blocks?  Not really...
hold on
plot(nanmedian(indPoint1), nanmean(indPoint2), 'o', 'markerSize', 14, 'markerFaceColor', 'r', 'markerEdgeColor', 'k', 'lineWidth', 1);
plot(indPoint1, indPoint2, '.')
plot([-10 2], [-10 2], '--k')
ylabel('indPoint2')
xlabel('indPoint1')

% summary of learning effects:
% there do not seem to be any effects of time on bias. Not by block, not by
% game. There are differences in draw fraction. games 2 and 3 seem to have
% lower drawing fraction than game 1. Game 2 makes sense... but game 3? Could be 
% because of timelimit but also could be learned. Learning interpretation
% is weakened by lack of block learning effects in any game.

% summary of value effects:
% there are not main effects of value on draw frequency or bias in games 1
% and 3. There are effects of value on draw frequency in game 2, as there
% should be. There also may be an interaction between value and bias
% whereby biased subjects become more biased at high values in games 2 and
% 3? Not clear, not enough data to tell whether this is real. 




% CHANGES ACROSS GAMES? 
% game 1
% stopSum1=stopSum
% infoFitBias1=infoFitBias;
% totDrawFrac1=totDrawFrac
% game 2
% stopSum2=stopSum;
% infoFitBias2=infoFitBias;
% totDrawFrac2=totDrawFrac;
% game 3
% stopSum3=stopSum;
% infoFitBias3=infoFitBias;
% totDrawFrac3=totDrawFrac;


% game 1 block 1
% stopSum4=stopSum
% infoFitBias4=infoFitBias;
% totDrawFrac4=totDrawFrac
% game 1 block 2
% stopSum5=stopSum
% infoFitBias5=infoFitBias;
% totDrawFrac5=totDrawFrac






% storing behavior of models with different inverse temperature settings
% modInfoFitBias1=modInfoFitBias;% invTemp=1
% modStopSum1=modStopSum;    % % invTemp=1      
% modInfoFitBias2=modInfoFitBias;% invTemp=10
% modStopSum2=modStopSum;    % % invTemp=10      
% 
% modInfoFitBias3=modInfoFitBias;% invTemp=100
% modStopSum3=modStopSum;    % % invTemp=100      
% 
% modInfoFitBias4=modInfoFitBias;% invTemp=.1
% modStopSum4=modStopSum;    % % invTemp=.1      
% 
% 
% modInfoFitBias5=modInfoFitBias;% invTemp=.5
% modStopSum5=modStopSum;    % % invTemp=.5    

% ideal:
% modStopSum = -7
% modInfoFitBias = -3.0187

x=modStopSum1;
y=modStopSum3;
minVal=min([x y])
maxVal=max([x y])
hold on
plot([minVal maxVal], [minVal maxVal], '--k')
plot(x, y, '.b')
plot(nanmean(x), nanmean(y), '*r');




