%% the goal of this script is to create some fake data, and
% then show that we can fit a model to recapture that data. 
% The model is a gaussian "information value" function that
% has 3 parameter:  mean (mu), stamdard deviation (sigma) and gain (gain). 
% the model chooses info based on a softmax function whereby value
% differences probabilisticaly translated into actions (ie bet or draw). 

% the hope is that this fitting procedure can be used on our behavioral
% data to assess the mean bias in information value for each subject.  We
% could also get other parameters... but this one seems like the first
% thing to focus on. 

% list of required file:
% infoValueFitTestCode.m
% softMax.m
% fitInfoValue.m




%% OK, step by step. 

% first, lets create a dataset like those in our task:

% starting with bead differences
xData=[-10:2:10 -10:2:10 -10:2:10 -10:2:10];

% Then we will assume that our "simulated subject" chooses based on a
% gaussian information value function and a fixed (known) information cost.
% in addition, our subject will have some randomness in his actual choices,
% but that is modeled later on.

infCost = .1;
mu      = -5; %mean
sigma   = 4; %variance
gain    = 3; %gain
invTemp = 100; %inverse temperature

% ok, now we can compute the "TRUE" inforamtion value function for this
% simulated subject:
infVal  = (normpdf(xData, mu, sigma).*gain)-infCost;

% reps is the number of times we will simulate data from the model
% and nChecks is the number of times we'll run the fitting procedure (with
% random start values, and choose the best ones).
reps=20;
nChecks=10;
for k = 1:reps
    
    %% ok, now lets have the subject do the task... he'll need to make some 
    %choices (here we'll just model inf/no info choices)
    pick=nan(length(infVal),1);
    for i = 1:length(infVal)
        pAll=softmax([infVal(i) 0], invTemp);
        pick(i) = find(cumsum(pAll)>rand, 1);
        % so pick==1 means chose info, pick==2 means no info
    end
    
    
    %% now lets see if we can fit a model to our simulated subject to recapture the parameters
    % that we used to create it. 
    
    % in order to make sure that the gradient descent algorithm finds the
    % right minimum, we're going to run through this a few times, but only
    % keep the parameters from the iteration where the model fit best
    for tt=1:nChecks
        infChoice=pick==1; % send in info choice trials
        % startpoint chooses some parameter values for the search to start
        % at
        startPoint= [(rand.*20)-10 rand.*3+1 1+rand.*5 30+rand.*10];
        % then we run a fitting function to choose the parameter that
        % minimize the negative log likelihood of parameters
        [c_fits, c_sse] = fitInfoValue(xData, infChoice, startPoint);
        % if this model fits better than the other ones we've tried, then
        % store the parameters
        if tt==1 ||c_sse<sse
            fits=c_fits;
            sse= c_sse;
        end
        % otherwise, try to fit it again with new random start conditions.
    end
    

    % make descriptive plot that we've been using to look at data
    possX=unique(xData);
    for i =1:length(possX)
        sel=xData==possX(i);
        drawFrac(i)=nanmean(infChoice(sel));
    end
    
    % compare that plot to the value function best fit to it. 
    
    hold on
    plot(possX, drawFrac, 'g')
    X=-10:10;
    infY=(normpdf(X, fits(1), fits(2)).*fits(3))-infCost;
    plot(X, infY./max(infY), 'r')
    infEst(k,:)=fits;
end

% if you ran the model a bunch of times, see how well the initial mean bias
% is recaptured:
close all
hist(infEst(:, 1))
