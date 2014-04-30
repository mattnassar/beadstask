%% OK, this script is just to calculate percent correct after a fixed
%% number of urn draws.
%% note to self... 2 beads at a time?

%% question:  do subjects tend to stop at an odd number of draws?
% should i make draws come 2 at a time?


q=.6

inR=0
inB=0
pRedIn=1./(1+(q./(1-q)).^((inR+inB)-(2.*(inR))));
pBlueIn=1-pRedIn;

inPctCorrect=max([pRedIn pBlueIn])


for i = 1:100

draws=i.*2;
totObs=inR+inB+draws;

%% Calculate probability distribution

newRedDraws=0:draws;
pRed=1./(1+(q./(1-q)).^((inR+inB+draws)-(2.*(inR+(newRedDraws)))));
pBlue=1-pRed;
weights = pRedIn.*binopdf(0:draws,draws, q)+pBlueIn.*binopdf(0:draws,draws, 1-q)
beadDifferential=    newRedDraws+inR-(totObs-newRedDraws-inR);
pctCorrect=max([pRed; 1-pRed])


% hold on
% plot(beadDifferential, weights, '--k')
% plot(beadDifferential, pctCorrect, 'k', 'lineWidth',10)
% plot(beadDifferential, pRed, 'r')
% plot(beadDifferential, 1-pRed, 'b')
% ylabel('pct correct')
% xlabel('bead differential')
% close all
 
expPctCorrect(i)=sum(pctCorrect.* weights)

end
