function [drawprob] = fitcurve(alpha, b, cv, invTemp, HP, LP)

% following code was copied and modified from statevaluecalc.m
limit = 30;
inc = 2;
probabilityprecision = 0;

depthlimit = 10;
calclowbound = -30;
calchighbound = 30;

UP = .60;
IC = 0.1;
HP = HP;
LP = LP;
HIP = 0;
LIP = 0;

x=1;

func1 = nan(1,2*limit/inc+1);
func2 = nan(1,2*limit/inc+1);
func3 = nan(1,2*limit/inc+1);
for BD = -limit:inc:limit
    [func1(x), func2(x), func3(x)] = calculator(BD, UP, IC, HP, LP, HIP, LIP, limit, 1, alpha, b, cv, probabilityprecision, depthlimit, 0, calclowbound, calchighbound);
    % bet high, bet low, draw
    x/length(-limit:inc:limit);
    x=x+1;
end
truncdraw = func3(limit/inc-5: limit/inc+7);
truncbethigh = func1(limit/inc-5: limit/inc+7);
truncbetlow = func2(limit/inc-5: limit/inc+7);
if length(truncdraw) ~= length(truncbethigh) || length(truncdraw) ~= length(truncbetlow)
    disp('error1')
    keyboard
end
drawprob = nan(1,length(truncdraw));
for eminem = 1:length(truncdraw)
    % softmax between the three
    drawprob(eminem) = 1/(1+exp((truncbethigh(eminem)-truncdraw(eminem))/invTemp)+exp((truncbetlow(eminem)-truncdraw(eminem))/invTemp));
end
end