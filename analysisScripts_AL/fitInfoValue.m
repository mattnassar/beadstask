function [estimates, sse] = fitInfoValue(xdata, ydata, start_point)

% x data should be bead differential (positive values favoring high value
% option)...

% y data should be binary indicating info choices

% start_point indicates parameters from which search begins


%% Goal:  fit a gaussian information value function and uniform draw cost function
%         to the beads task data. 
%  parameters:  mean, variance, gain, softmax inverse temperature
%


if nargin<3
start_point = [0, 1, 10, 100]; % this supplies the initial parameters for the curve search
end

model = @infFun;                % this is the function handle to the function that takes the parameters and outputs the thing we want to minimize

options=optimset('MaxFunEvals', 1000000);
[estimates, sse, ef, o, l, g, h] = fmincon(model, start_point, [], [], [], [], [-20,  1, 0 , 0],[20, 10e20, 10e20, 10e20], [], options );

 
    function [negLogLike, FittedCurve] = infFun(params)  %input comes from fmincon (initially start_point)
        infCost=.1;  % inf cost is fixed for now...
       % keyboard
        mu      = params(1); %mean
        sigma   = params(2); %variance
        gain    = params(3); %gain
        invTemp = params(4); %inverse temperature
        infVal  = (normpdf(xdata, mu, sigma).*gain)-infCost;
        
            for i = 1:length(infVal);
               pAll(i,:)=softmax([infVal(i) 0], invTemp);
            end
        
          %  keyboard
            
       logLike=sum(log(pAll(ydata, 1)))+ sum(log(pAll(~ydata, 2)));
       negLogLike = logLike.*-1   ;
       
       if ~isfinite(negLogLike)
           keyboard
       end
       
       
       
    end
end
