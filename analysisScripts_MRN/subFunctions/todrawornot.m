function [draw] = todrawornot(beaddiff)

% This function takes beaddiff and spits out some drawing behavior based on
% previous subject choices. 


if beaddiff >= 8 | beaddiff <= -12
    draw = 0;
else
%     measured=  [0 .1786 0.3103 0.4194 0.5517 0.4138 0.1724 0.1613 0.1 0.0938 0]
   % x       = -12:2:12
    interpX=  -12:1:12;
%     interpF=  interp(measured, 2)
%     interpF(interpF<0) = 0;
    % hard code to make it quicker.
    F=[ 0    0.0935    0.1786    0.2523    0.3103    0.3572    0.4194 ...   
    0.5017    0.5517    0.5211    0.4138    0.2735    0.1724    0.1528 ...
      0.1613    0.1346    0.1000    0.0962    0.0938    0.0587   0        0];    
    
    frac=F(interpX==beaddiff);

    %rng('shuffle')
    if rand(1) <= frac
        draw = 1;
    else
        draw = 0;
    end
end

