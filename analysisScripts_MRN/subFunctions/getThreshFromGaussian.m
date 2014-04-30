function [threshEsts]=getThreshFromGaussian(gaussParams)      

% basically we are just going to solve for the points at which this curve
% crosses y=1./3
 
U=    gaussParams(1); %mean
ss=    gaussParams(2); %standard deviation
gain=  gaussParams(3); %gain
offset =gaussParams(4);% offset

if offset>(1./3)
    threshEsts=[-inf inf];
else




% we will solve the equation where the pdf of a gaussian =1./3 for x.
% this equation boils down to a quadratic where A=1, B=-2U,
% and...
C=(2.*ss.^2).*log((sqrt(2.*pi.*ss).*((1./3)-offset))./(gain))-U.^2;
A=1;
B=-2.*U;

plusHalf=(-B+sqrt((B.^2)-(4.*A.*C)))./(2.*A);
minusHalf=(-B-sqrt((B.^2)-(4.*A.*C)))./(2.*A);
threshEsts=[plusHalf minusHalf];

end




