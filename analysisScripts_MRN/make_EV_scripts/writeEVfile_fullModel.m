function [regNames timings]=writeEVfile_fullModel(timingData, allRegModVars, evFileName);



%% OK, what is the goal here:

% here we'll create the EV files for a basic model based analysis


% THE BASIC ANALYSIS:
% model timings:

% trial viewing 
% info viewing
% button push (draw, noDraw, bet(left or right) all modeled together)
% bet (lock in bet with button push)
% "too slow" -- tooSlowOnTime + 1 second
% "force bet" -- forceChoiceTime + .5 seconds

% One tricky thing will bet getting the times of "pushed bet".
% it should be enterDecTime, but only on trials where forceChoiceTime is
% nan. Otherwise there was no button pushed, as the bet was "forced". 

% model modulators:

% trial viewing -- 1) number of beads, 2) bead difference?
% info viewing  -- 1) info value, 2) info gain, 3) KL divergence
% bet -- 1) was it left? 

% extra button push -- 1) was it left?, 2) was it extra draw or betNow?
% 


% The following analysis will come second, and we can check the validity of
% our trial coefficients by comparing results to the first model. I guess
% we can also check validity by testing face/house responses in PPA/FFA.

% THE MONEY ANALSYSIS:
% the money analysis will require estimating trial information
% coefficients. Once we have these coefficients we can use them to 1) try
% to classify houses/faces... then test how trial classification scores 
% depend on information value, information gain (decreased entropy), and KL
% divergence. 




%keyboard



%% make beads task EV files.  MRN 2-12-14
% send in timing and model based data for a run... this function will write
% an EV file for each regressor and save it as "evFileName".

% timingData already has the scan onset time subtracted, meaning it should
% be perfectly aligned for FSL.

% however, we may want to drop a few frames at the beginning of the scan.
dropFrames=5;
TR=1.5;
timeDropped=dropFrames.*TR;
ll=length(allRegModVars.info.vals);
% timing variables require an "on" timestamp and a duration. In order to
% use modulator regressors, each timing variable must be matched in name to
% the fieldnames of the allRegModVars structure.


% apparently the last allRegModVars.info column used to be the house/face
% variable. that is currently not the case though...

% % the goal of this model is to look at info related variables directly, and
% % since these tend to correlate within block with the house/face variable,
% % i'll not include that variable in this model.
% allRegModVars.info.vals=(allRegModVars.info.vals(:,1:end-1))
% allRegModVars.info.names=allRegModVars.info.names(1:end-1);
% 
%keyboard

timings(1).on          =timingData.choiceOn;
timings(1).duration    =nanmedian(timingData.choiceOff-timingData.choiceOn);
timings(1).name        ='choice';

timings(2).on          =timingData.infoOn;
timings(2).duration    =nanmedian(timingData.infoOff-timingData.infoOn);
timings(2).name        ='info';

timings(3).on           =timingData.betOn;
timings(3).duration     =nanmedian  (timingData.betChoiceTime- timingData.betOn);
timings(3).name        ='bet';

timings(4).on           =timingData.forceChoiceTime;
timings(4).duration     =.5;
timings(4).name        ='forceBet';


% OK, there seems to be an issue with these... lets hold off until arthur
% looks into the timestamps a bit...

% timings(5).on           =timingData.tooSlow;
% timings(5).duration     =1;
% timings(5).name         ='tooSlow';



% model extra button pushes. 
% we need: 1) button push timings, 2) info or bet push, 3) right or left
% push. 

timings(5).duration     =.1;
timings(5).name         ='buttonPush';
timings(5).on           =[]

isDraw=[];
allPushSide=[];
for i = 1:ll
     if ~isfinite(timingData.forceChoiceTime(i))
        buttonPushes=[timingData.extraButtonPush{i}, timingData.betOn2(i)];
        trialDraw=[true(length(timingData.extraButtonPush{i}),1); false];
     else
        buttonPushes=[timingData.extraButtonPush{i}];
        trialDraw=[true(length(timingData.extraButtonPush{i}),1)];
     end
     % was the right or left button pushed? decode from info button and
     % info push
     pushSide = (trialDraw.*2 - 1).*(allRegModVars.infOnRight(i).*2 - 1);

    buttonPushes=buttonPushes(isfinite(buttonPushes));
    trialDraw   =trialDraw(isfinite(buttonPushes));
    pushSide    =pushSide(isfinite(buttonPushes));
    
    isDraw=[isDraw; trialDraw;];    
    allPushSide=[allPushSide; pushSide]; 
    
    trialExPush=buttonPushes;
    timings(5).on=[ timings(5).on; trialExPush'];
end
%timings(6).on=timings(6).on(isfinite(timings(6).on));

allRegModVars.buttonPush.names={'draw', 'right'};
allRegModVars.buttonPush.vals =meanCentX([isDraw, allPushSide]);
rmfield(allRegModVars, 'infOnRight') 




k=1;

modTimes=fieldnames(allRegModVars);
clear regNames
for i = 1:length(timings)
    
    sel=isfinite(timings(i).on);
    
    % create 3 column matrix
    featTxt=[timings(i).on(sel)-timeDropped, ...                    % onset
        ones(sum(sel),1).*timings(i).duration, ... % duration
        ones(sum(sel),1)];                         % height
    
    if isempty(featTxt)
        featTxt=[0 1 0]; % send in some bullshit if its empty.  Maybe fsl will know what to do...
    end


    % write text file
    dlmwrite([evFileName  '_reg_' num2str(k)], featTxt, ' ')
    regNames{k}=[timings(i).name '_timing']; % store regressor name
    k=k+1; % increment regressor number.
    
    % if there are modulators, we'll want to make EV files for them too.
    if ~isempty(strmatch(timings(i).name, fieldnames(allRegModVars)))

        % get the modulators of this particular epoch.
        eval(sprintf('timeMods=allRegModVars.%s', timings(i).name))
        
        % check to make sure that timing and modulator arrays are same
        % length:
         
        if length(timeMods.vals)~=length(timings(i).on)
            disp('houston, we have a problem');
            keyboard;
        end
        %keyboard
        % run through each modulator of this epoch and create a file.
        for j = 1:length(timeMods.names)
            featTxt=[timings(i).on(sel)-timeDropped, ...              % onset
                ones(sum(sel),1).*timings(i).duration, ...  %duration                  % duration
                ones(sum(sel),1).*timeMods.vals(sel,j)];                   % height
            
            dlmwrite([evFileName  '_reg_' num2str(k)], featTxt, ' ')
            regNames{k}=[timings(i).name timeMods.names{j}]; % store regressor name
            k=k+1;
        end
    end
end




