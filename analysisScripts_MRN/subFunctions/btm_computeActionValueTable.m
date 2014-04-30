function [drawValue pickRedVal pickBlueVal blueThresh redThresh pBlueUrn, pBlueBallOnNextDraw] ...
    = btm_computeActionValueTable(treeHeight, q, infCost, bCorr, rCorr, bErr, rErr)


%% bead task model (STRIPPED!!!)


% treeHeight=200;
% q=.6;    % this is the probability of the more likely ball draw
% infCost=.1;
% 
% bCorr=100;
% rCorr=10;
% bErr=0;
% rErr=0;



%% create tree
b=0;  % this is to manipulate start position (blue)
r=0 ; % same deal... but red

%% compute initial parameters:
pBgiv_br=1./(1+(q./(1-q)).^((b+r)-(2.*b)));
pRgiv_br=1-pBgiv_br;
%pBball=pBgiv_br.*q+pRgiv_br.*(1-q);

%% compute probability of getting to all states in the tree.
%% OK, 5-3-13... lets do the tree a little different. 
possReds=repmat(0:treeHeight, treeHeight+1, 1);
possBlues=repmat((0:treeHeight)', 1, treeHeight+1);
impossible=possReds+possBlues>treeHeight;
possReds(impossible)=nan;
possBlues(impossible)=nan;

%% compute the probability of the urn being blue at each of those states:
pBlueUrn=1./(1+(q./(1-q)).^((possBlues+possReds+b+r)-(2.*(possBlues+b))));
pRedUrn=1-pBlueUrn;


%% compute the probability of drawing bead-types from each state:

pBlueBallOnNextDraw=(pBlueUrn.*q)+pRedUrn.*(1-q);
pRedBallOnNextDraw=(pBlueUrn.*(1-q))+pRedUrn.*(q);


%% get urn picking values

% this is the expected value of taking the action: Draw.
pickRedVal=pRedUrn.*rCorr + (1-pRedUrn).*rErr;
pickBlueVal=pBlueUrn.*bCorr + (1-pBlueUrn).*bErr;
pickVal=  max(cat(3, pickBlueVal, pickRedVal), [], 3);

% this is the expected value of the state that you'll be in if you drew a
% bead.
drawValue=nans(size(pBlueBallOnNextDraw));
% this is the maximum of the two possible values.
stateValue=nans(size(pBlueBallOnNextDraw));

for i=treeHeight+1:-1:1
    sel=possReds+possBlues==i-1;
    if i==treeHeight+1
        drawValue(sel)=-inf;
    else
        [expStateVal]=computeFutureStateExpectation(stateValue, pBlueBallOnNextDraw, pRedBallOnNextDraw);
        smSel=sel(1:end-1,1:end-1);
        drawValue(sel)=expStateVal(smSel)-infCost;
    end
    stateValue=max(cat(3, drawValue, pickVal), [], 3);
end


actionDiff=drawValue-pickVal;
find((actionDiff>0));
[I,J] = find(actionDiff>0);
blueThresh=min(J-I).*-1;
redThresh=max(J-I);
