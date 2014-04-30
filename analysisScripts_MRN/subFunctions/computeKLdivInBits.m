function [klDiv]=computeKLdivInBits(discDist1, discDist2)
    % computes Dkl(P||Q) where discrete probability distribution P is a
    % vector (discDist1) and Q is another discrete probability distribution
    % on the same space (discDist2).
    
    % this notation is a bit funny, because for Bayesian updating you want
    % the posterior to be discDist1 and the prior to be discDist2... that
    % way you are computing how many extra bits you needed to encode the
    % data without the most recent update.
         
    
        
         p=discDist1.*log2(discDist2./discDist1);
         p(discDist1==0)=0;
         klDiv=-1.*nansum(p);
         
         
         
         
         
         
         
         
         
         