function [discStat, inSampAcc, pVal, unbiasedAccMeasure, subB]=getDiscStat(yVar, xMat, numPerms)
    

%    THIS CODE WAS FOR CHECKING GLM FIT.  
%  it was behaving suspiciously, but it seems to check out, at least for
%  large datasets...
%     xMat=normrnd(0, 1, 1000000, 2);
%     yProb=1./(1+ exp(-xMat(:,1)));
%     yVar=binornd(1,yProb);

if nargin<3||isempty(numPerms)
    numPerms=100;
end



[subB, dev]=glmfit([xMat], yVar,'binomial');
testStat=(subB'*[ones(length(xMat),1), xMat]')';
discStat=testStat;
discStat(~logical(yVar))=testStat(~logical(yVar)).*-1;
inSampAcc=nanmean(discStat>0);

if numPerms>0
    for i = 1:numPerms
        permYs=yVar(randperm(length(yVar)));
        [B, dev]=glmfit([xMat], permYs,'binomial');
        ptestStat=(B'*[ones(length(xMat),1), xMat]')';
        pdiscStat=ptestStat;
        pdiscStat(~logical(yVar))=testStat(~logical(yVar)).*-1;
        pInSampAcc(i)=nanmean(pdiscStat>0);
    end
    pVal=sum(inSampAcc>pInSampAcc)./numPerms;
end
unbiasedAccMeasure=inSampAcc-nanmedian(pInSampAcc);

