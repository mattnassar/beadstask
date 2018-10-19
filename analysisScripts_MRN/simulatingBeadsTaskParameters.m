
%

addpath(genpath('.'))

%
%
% - jar ratio of blue:orange beads: 60/40, 65/35, 55/45
% - payoffs: 70/10, 80/10, 60/10
% - cost to draw beads: .10, .15, .20, .50
% - starting beads: -10:10, -5:5
%
%


varStruct.majBeadProp=[.65, .60, .55];
varStruct.maxPayout=[70, 80, ];

% Loop through a bunch of starting conditions and see what optimal behavior
% should look like:
maxPayouts =[ 70]
drawCosts  =[.1];
majBeadPros=[.6];
startBeads=-10:10;




for ii = 1:length(maxPayouts)
    for jj = 1:length(drawCosts)
        for kk = 1:length(majBeadPros)
            
            
            maxPayout =maxPayouts(ii);
            drawCost  =drawCosts(jj);
            majBeadPro=majBeadPros(kk);
            
            
            
            lowValPayout=10;
            
            
            
            % trialData
            % maxUrnProb
            % startLeft
            % startRight
            % hiValSide (0 = left, 1 = right)
            % infCost
            % hiValue
            % loValue
            % inValue
            
            close all
            
            allPenalties=-100:0;
            for i=1:length(allPenalties)
                
                trialData.startLeft =0:40;
                trialData.startRight=ones(size(trialData.startLeft)).*20;
                trialData.maxUrnProb =ones(size(trialData.startLeft)).*majBeadPro;
                trialData.hiValSide =true(size(trialData.startLeft));
                trialData.infCost =ones(size(trialData.startLeft)).*drawCost;
                trialData.hiValue =ones(size(trialData.startLeft)).*maxPayout;
                trialData.loValue =ones(size(trialData.startLeft)).*lowValPayout;
                trialData.inValue =ones(size(trialData.startLeft)).*allPenalties(i);
                trialData.corrUrnSide=true(size(trialData.startLeft));
                trialData=straightStruct(trialData);
                [infoData, modBehav, condInfoMetrics]=getInfoValueFromBlock(trialData, 100)
                
                doDraw(:,i)=(modBehav.choiceProbs(:,3)>.5);
                
                
            end
            
            evidence=trialData.startRight-trialData.startLeft;
            
            figure(1)
            
            
            hold on
            imagesc(allPenalties, evidence,  doDraw, [0, 1])
            plot([-100, 0], [0, 0], '--k')
            ylabel('evidence for high value side')
            xlabel('penalty')
            
            
            fn=sprintf('penaltyEffects_%g_%g_%g.eps', maxPayout, majBeadPro,drawCost );
            ylim([-10, 10])
            saveas(gcf, fn,'epsc2')
            close(gcf)
            
            
            
            figure(2)
            
            hold on
            ylabel('optimal choice prob')
            xlabel('evidence for high value option')
            
            plot(trialData.startRight-trialData.startLeft,   infoData.modTrialInfValue./max(infoData.modTrialInfValue), 'g')
            ff=legend('choose low', 'choose high', 'draw', 'information value');
            plot([0, 0], [0, 1], '--k')
                        plot([-20, 20], [0, 0], '--k')

            set(ff, 'box','off')
            set(gca, 'box','off')
            
            xlim([-10, 10])
            fn=sprintf('modBehav_penalty_%g_%g_%g.eps', maxPayout, majBeadPro,drawCost );
            
            saveas(gcf, fn, 'epsc2')
            
            close(gcf)
            
        end
    end
end



