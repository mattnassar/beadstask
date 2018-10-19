%% load all scanner data
%% get everything we need on the path.
clear classes

rootDir='/Users/mattnassar/matt_work_stuff/Matt/m_files/beadsTask/analysisScripts_MRN/'
cd(rootDir);
addpath(genpath(rootDir))


% TO DO:
% 1) get some measures of choice bias 
% 2) look at choice bias versus draw bias
% 3) why are some people not drawing? Look block by block?


%% load subject data
%subFileDir='/Users/mattnassar/Dropbox/BeadsTaskCode/BeadsTask4';

subFileDir='/Users/mattnassar/Dropbox/beadstaskcode/BeadsTask5data';
figDir='/Users/mattnassar/Dropbox/BeadsTaskCode/figures/'
subNames={'CC528', 'FL1136', 'FT3594', 'IA3593', 'IL3520',	'LL3555', ...
    'NS950', 'OF3592', 'TK3556', 'TQ3543', 'TQ3600', ...
    'BA3550', 'KQ3548', 'LL3555', 'MF3567',  'NN3554' }


% DATA FROM OUT OF SCANNER PILOT FOR FORCED DRAW VERSION:
% subNames={'BE1829',	'KK1101',	'testbug'};


% DATA FROM INITIAL fMRI PILOT:
% subNames={'DZ572',	'HE3515',	'KQ3513',	'LT3516',	'LZ3512',	'NN3511', ...
%   'OK547',	'SI3517',   'EO3518', 	'HU2138',	'LL2026',	'LX1803', ...
%   'MC1906',	'NS969',	'QB2027',	'TA3514'};

% loop through each subject and load data, store in "allSubjData" structure
% and save that structure so we don't need to use this slow loading code
% again.

for i = 1:length(subNames)
    fileName=fullfile(subFileDir, subNames{i});

    
    [allData]=beadTaskDataLoader_mrn(fileName);

    
    % put data structures in a bigger structure including all subjects.
    eval(sprintf('allSubjData.%s=allData;', subNames{i}));
    
    % beadTaskDataLoader changes directories, so lets go back manually.
    cd(subFileDir)
end


cd(figDir)
% 'game1.block1', 'game1.block2',   'game1', 'game1',

noDrawBlkNames={'noDraw_asym.block1', 'noDraw_asym.block2'};
blkNames={ 'realScannerGame1.block1', 'realScannerGame1.block2', 'realScannerGame2.block1','realScannerGame2.block2'}; %
blk={ 'realScannerGame1', 'realScannerGame1', 'realScannerGame2', 'realScannerGame2'};
behavSubjs=subNames;





% 
% % Testing unpackBeadsData:
% [allSubjData.IL3520.realScannerGame1.block1.statusData.extraDrawTime]
% 
% 
% data=unpackBeadsData(allSubjData.IL3520.realScannerGame1.block1.statusData)



makeFigs=false;
invT=1;
ll=length(subNames)
biasSum=nan(ll,1);
biasMax=nan(ll,1);
stopSum=nan(ll,1);
hiValStopSum  =nan(ll,1);
loValStopSum  =nan(ll,1);
hiValindPoint =nan(ll,1);
loValindPoint =nan(ll,1);
hiValDrawFrac =nan(ll,1);
loValDrawFrac =nan(ll,1);
totDrawFrac   =nan(ll,1);
totModDrawFrac=nan(ll,1);
totThreshDist =nan(ll,1);
totModThreshDist=nan(ll,1);
noDrawBias=nan(ll,1);
indPoint=nan(ll,1);
infoFitBias=nan(ll,1);
stopFitSum=nan(ll,1);
stopMedSum=nan(ll,1);

nonUndergradVar=nan(ll,1);
infoFitStd=nan(ll,1);
avThresh=nan(ll,1);
modInfoFitBias=nan(ll,1);
modStopSum=nan(ll,1);
defaultPlotParameters
drawThresh=1;

for i = 1:length(behavSubjs)    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              real scanner task analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                exValues=straightStruct(exValues);
                
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
        nDraws=nan(size(allData.curr_trialInfList));
        for k=1:length(allData.curr_trialInfList);
            nDraws(k)=length(allData.curr_trialInfList{k})
        end
        subDraw=nDraws>drawThresh;
         
        
        diffs=unique(allData.begBeadDiff);
        for k = 1:length(diffs);
            sel=allData.begBeadDiff==diffs(k);
            meanDrawFrac(k)=nanmean(subDraw(sel));
            meanBetHiFrac(k)=nanmean(~subDraw(sel)&(allData.hiValSide(sel)==allData.curr_choice(sel)));
            meanBetLoFrac(k)=nanmean(~subDraw(sel)&(allData.hiValSide(sel)~=allData.curr_choice(sel)));
            
            meanModDrawFrac(k)=nanmean(allModData.choiceCategory(sel)==3);
            meanModBetHiFrac(k)=nanmean(allModData.choiceCategory(sel)==2);
            meanModBetLoFrac(k)=nanmean(allModData.choiceCategory(sel)==1);
        end
          

        
        if makeFigs
            close all
            hold on
            plot([indPoint(i) indPoint(i)], [0 1], 'c')
            plot(diffs, meanDrawFrac, 'g')
            plot(diffs, meanBetHiFrac, 'r')
            plot(diffs, meanBetLoFrac, 'b')
            plot([0 0], [0 1], '--k')
            xlabel('Bead difference')
            ylabel('Choice frequency')
            %plot(diffs, fitProb, '--m')
            
            set(gca, 'box', 'off')
            saveas(gcf, sprintf('indSubDrawPlot_%s.eps', behavSubjs{i}), 'epsc2')
            %keyboard
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
        %plot(diffs, modFitProb, '--m')
        xlabel('Bead difference')
        ylabel('Choice frequency')
        set(gca, 'box', 'off')
        close all
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  map thresholds based on end of betting: 
        %  HERE WE WANT TO DROP TRIALS WHERE THE SUBJECT DID NOT
        %  INTENTIONALLY DRAW A BEAD. THIS MEANS THAT NDRAWS SHOULD BE
        %  GREATER THAN 1.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
         
        allPossDiffs=-30:30;
        sel=subDraw&(allData.hiValSide==allData.curr_choice);
        if sum(sel)>1
            hCount=histc(allData.endBeadDiff(sel), allPossDiffs);
            [estimates, sse] = fitGaussPDF_wGain (allPossDiffs', hCount, []);
            hDistFit=estimates(1);
            hDistMed=median(allData.endBeadDiff(sel));
        else
            hDistFit=nan;
            hDistMed=nan;
        end
        
        
        sel2=subDraw&(allData.hiValSide~=allData.curr_choice);
        if sum(sel2)>1
            lCount=histc(allData.endBeadDiff(sel2), allPossDiffs);
            [estimates, sse] = fitGaussPDF_wGain (allPossDiffs', lCount, []);
            lDistFit=estimates(1);
            lDistMed=median(allData.endBeadDiff(sel2));

        else
            lDistFit=nan;
            lDistMed=nan;
        end
        
        stopFitSum(i)=hDistFit+lDistFit;
        stopMedSum(i)=hDistMed+lDistMed;        
        highStop = nanmean(allData.endBeadDiff(subDraw&(allData.hiValSide==allData.curr_choice)));
        lowStop  = nanmean(allData.endBeadDiff(subDraw&(allData.hiValSide~=allData.curr_choice)));
        % compute same thing for high and low value trials separately
        hiHighStop = nanmean(allData.endBeadDiff(isHi&subDraw&(allData.hiValSide==allData.curr_choice)));
        loHighStop = nanmean(allData.endBeadDiff(~isHi&subDraw&(allData.hiValSide==allData.curr_choice)));
        hiLowStop  = nanmean(allData.endBeadDiff(isHi&subDraw&(allData.hiValSide~=allData.curr_choice)));
        loLowStop  = nanmean(allData.endBeadDiff(~isHi&subDraw&(allData.hiValSide~=allData.curr_choice)));
        
        modHighStop= nanmean(allModData.endDiff(allModData.choiceCategory==3&(allModData.betSide==2)));
        modLowStop= nanmean(allModData.endDiff(allModData.choiceCategory==3&(allModData.betSide==1)));
        



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

mx=max(abs([stopSum; modStopSum])) +.25;
my=max([avThresh; totModThreshDist])+.25;
minY=min([avThresh; totModThreshDist; 0])-.25;

hold on
plot([0 0], [0 my], '--k')
a=plot(stopSum, avThresh, 'o', 'markerSize', 10, 'markerFaceColor', cbColors(2,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
b=plot(modStopSum, totModThreshDist, 'o', 'markerSize', 10, 'markerFaceColor', cbColors(3,:), 'markerEdgeColor', 'k', 'lineWidth', 1);
aa=legend([a, b], 'subjects', 'model')
set(aa, 'box', 'off')
set(gca, 'box', 'off')
xlim([-1, 1].*mx);
ylim([minY, my]);
ylabel('average threshold')
xlabel('threshold bias')
saveas(gcf, 'threshVsBiasPlot.eps', 'epsc2')
close all


