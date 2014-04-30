function [estimates, sse] = fitGaussPDF_wGain (xdata, ydata, startPoint)
%% fitGaussPDF_wGain, MRN, 1/20/2014.
%  created fom shell of fitGaussPDF, MRN, who knows.

%  this code should be capable of fitting a guassian distribution  across
%  input (xdata) to output (ydata). Here we will assume that input 



if nargin<3|isempty(startPoint)
start_point = [nansum(xdata.* ydata./nansum(ydata)), ...  MEAN
              nanstd(xdata), ...   STD
              nanmean(ydata), ... GAIN
               0];              %  OFFSET
end

model = @GaussFun;                % this is the function handle to the function that takes the parameters and outputs the thing we want to minimize

options=optimset('MaxFunEvals', 10e30);
[estimates, sse, ef, o, l, g, h] = fmincon(model, start_point, [], [], [], [], [-60, .00001, 0, 0],[60, 10000, 10000, 1], [], options );

 
    function [sse, FittedCurve] = GaussFun(params)  %input comes from fmincon (initially start_point)
        U =    params(1); %mean
        ss=    params(2); %standard deviation        
        gain=  params(3); %gain
        offset =params(4);% offset
        FittedCurve =offset+(gain.*(1/sqrt(2*pi*ss) * exp(-(xdata-U).^2/(2*ss))));
        FittedCurve(~isfinite(FittedCurve))=10e1000;        
        ErrorVector = FittedCurve - ydata  ;                    %error at each point along the curve
        sse = nansum(ErrorVector.^2) ;                        %SSE= output= what will be minimized
        
    end
end

