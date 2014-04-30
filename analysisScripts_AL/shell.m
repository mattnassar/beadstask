%% import data
%#ok<*UNRCH>  % this creates exception for unreachable codes inside if-0 statements
clear classes
clc
cd /Users/Carmenere/Documents/MATLAB
% subjects for which you have behavioral data
list = [{'NN3511'};{'LZ3512'};{'KQ3513'};{'DZ572'};...
    {'TA3514'};{'HE3515'};{'QB2027'};{'LT3516'};...
    {'OK547'};{'LX1803'};{'MC1906'};{'NS969'};...
    {'SI3517'};{'HU2138'};{'LL2026'};{'EO3518'};...
    {'BJ2376'};{'BX2031'};{'DT2005'};{'NC3509'};{'SQ3510'}]; %%these five were pilots

% pilot subjects
% list = [{'BJ2376'};{'BX2031'};{'DT2005'};{'NC3509'};{'SQ3510'}];
sump = zeros(1,11);
sumt = zeros(1,11);
str = cd;
totalData = struct;
for m = 1:length(list)
    subjID = list{m,:};
    [allData] = beadTaskDataLoader(subjID);
    eval(sprintf('totalData.%s=allData', subjID));
end
cd(str)
save totalData.mat totalData

%% analyze drawing
if 1
    clear classes
    clc
    load totalData.mat
    allNames = fieldnames(totalData);
    choices = struct;
    for m = 1:length(allNames)
        lowreward = struct;
        lowreward.drewbeads = zeros(1,21);
        lowreward.pickedlow = zeros(1,21);
        lowreward.pickedhigh = zeros(1,21);
        lowreward.correctbet = zeros(1,21);
        lowreward.incorrectbet = zeros(1,21);
        lowreward.total = zeros(1,21);
        highreward = struct;
        highreward.drewbeads = zeros(1,21);
        highreward.pickedlow = zeros(1,21);
        highreward.pickedhigh = zeros(1,21);
        highreward.correctbet = zeros(1,21);
        highreward.incorrectbet = zeros(1,21);
        highreward.total = zeros(1,21);
        for n = 1:2
            if n == 1
                statusData = totalData.(allNames{m}).scan.block1.statusData;
            elseif n == 2
                statusData = totalData.(allNames{m}).scan.block2.statusData;
            else
                disp('Error Code 1');
                keyboard
            end
            
            for q= 1:length(statusData)
                leftishigh = 0;
                rightishigh = 0;
                if statusData(q).isGoodTrial == 1
                    if statusData(q).curr_rewCorrLeft > statusData(q).curr_rewCorrRight
                        % left option has higher payout
                        leftishigh = 1;
                        temp = statusData(q).curr_start_tokensLeft - statusData(q).curr_start_tokensRight;
                        % temp is the number of beads favoring the high payout option
                    elseif statusData(q).curr_rewCorrLeft < statusData(q).curr_rewCorrRight
                        % right option has higer payout
                        rightishigh = 1;
                        temp = statusData(q).curr_start_tokensRight - statusData(q).curr_start_tokensLeft;
                        % temp is the number of beads favoring the high payout option
                    elseif statusData(q).curr_rewCorrLeft == statusData(q).curr_rewCorrRight
                        % equivalent payout
                        rightishigh = 1;
                        temp = statusData(q).curr_start_tokensRight - statusData(q).curr_start_tokensLeft;
                    elseif isempty(statusData(q).curr_rewCorrLeft) || isempty(statusData(q).curr_rewCorrRight)
                        % empty cells
                    else
                        %unexpected error
                        disp('Error Code 3');
                        keyboard
                    end
                    
                    if statusData(q).curr_penErrLeft ~= 0 %|| (statusData(q).curr_rewCorrLeft + statusData(q).curr_rewCorrLeft >70)
                        highreward.total(temp+11) = highreward.total(temp+11)+1;
                        if ~isempty(statusData(q).curr_trialInfList)
                            highreward.drewbeads(temp+11) = highreward.drewbeads(temp+11)+1;
                        else
                            if (leftishigh && statusData(q).curr_choice == 1) || (rightishigh && statusData(q).curr_choice == 2)
                                highreward.pickedhigh(temp+11) = highreward.pickedhigh(temp+11)+1;
                            else
                                highreward.pickedlow(temp+11) = highreward.pickedlow(temp+11)+1;
                            end
                        end
                        if statusData(q).curr_choice == statusData(q).curr_urnType
                            highreward.correctbet(temp+11) = highreward.correctbet(temp+11)+1;
                        else
                            highreward.incorrectbet(temp+11) = highreward.incorrectbet(temp+11)+1;
                        end
                    else
                        lowreward.total(temp+11) = lowreward.total(temp+11)+1;
                        if ~isempty(statusData(q).curr_trialInfList)
                            lowreward.drewbeads(temp+11) = lowreward.drewbeads(temp+11)+1;
                        else
                            if (leftishigh && statusData(q).curr_choice == 1) || (rightishigh && statusData(q).curr_choice == 2)
                                lowreward.pickedhigh(temp+11) = lowreward.pickedhigh(temp+11)+1;
                            else
                                lowreward.pickedlow(temp+11) = lowreward.pickedlow(temp+11)+1;
                            end
                        end
                        if statusData(q).curr_choice == statusData(q).curr_urnType
                            lowreward.correctbet(temp+11) = lowreward.correctbet(temp+11)+1;
                        else
                            lowreward.incorrectbet(temp+11) = lowreward.incorrectbet(temp+11)+1;
                        end
                    end
                end
            end
        end
        lowreward.drewbeads = lowreward.drewbeads(1:2:21);
        lowreward.pickedlow = lowreward.pickedlow(1:2:21);
        lowreward.pickedhigh = lowreward.pickedhigh(1:2:21);
        lowreward.correctbet = lowreward.correctbet(1:2:21);
        lowreward.incorrectbet = lowreward.incorrectbet(1:2:21);
        lowreward.total = lowreward.total(1:2:21);
        highreward.drewbeads = highreward.drewbeads(1:2:21);
        highreward.pickedlow = highreward.pickedlow(1:2:21);
        highreward.pickedhigh = highreward.pickedhigh(1:2:21);
        highreward.correctbet = highreward.correctbet(1:2:21);
        highreward.incorrectbet = highreward.incorrectbet(1:2:21);
        highreward.total = highreward.total(1:2:21);
        % check
        if lowreward.drewbeads+lowreward.pickedlow+lowreward.pickedhigh ~= lowreward.total
            disp('check1 error');
            keyboard
        elseif highreward.drewbeads+highreward.pickedlow+highreward.pickedhigh ~= highreward.total
            disp('check2 error');
            keyboard
        elseif lowreward.correctbet+lowreward.incorrectbet ~= lowreward.total
            disp('check3 error');
            keyboard
        elseif highreward.correctbet + highreward.incorrectbet ~= highreward.total
            disp('check4 error');
            keyboard
        end
        eval(sprintf('%s=struct;', allNames{m}));
        eval(sprintf('%s.lowreward = lowreward;', allNames{m}));
        eval(sprintf('%s.highreward = highreward;', allNames{m}));
        eval(sprintf('choices.%s = %s;', allNames{m},allNames{m}));
    end
    save choices.mat choices
    disp('choices analyzed!');
end

%% plot some data
if 1
    clear classes
    load choices.mat
    undergrads = [{'NN3511'};{'LZ3512'};{'KQ3513'};{'DZ572'};{'TA3514'};{'HE3515'};...
        {'OK547'};{'NS969'};{'SI3517'};{'EO3518'};{'BJ2376'};{'NC3509'};{'SQ3510'}];
    nonundergrads = [{'LX1803'};{'QB2027'};{'HU2138'};{'LL2026'};{'MC1906'};{'LT3516'};{'BX2031'};{'DT2005'}];
    if 0
        for something = 1:length(nonundergrads)
            choices = rmfield(choices, nonundergrads{something});
        end
    end
    x = -10:2:10;
    allNames = fieldnames(choices);
    highdrewbeads = zeros(1,11);
    highpickedlow = zeros(1,11);
    highpickedhigh = zeros(1,11);
    highcorrect = zeros(1,11);
    lowdrewbeads = zeros(1,11);
    lowpickedlow = zeros(1,11);
    lowpickedhigh = zeros(1,11);
    lowcorrect = zeros(1,11);
    lowmean = nan(1,length(allNames));
    highmean = nan(1,length(allNames));
    lowsum = nan(1,length(allNames));
    highsum = nan(1,length(allNames));
    for q = 1:length(allNames)
        lownumdraws = sum(choices.(allNames{q}).lowreward.drewbeads);
        highnumdraws = sum(choices.(allNames{q}).highreward.drewbeads);
        choices.(allNames{q}).lowreward.drewbeads = choices.(allNames{q}).lowreward.drewbeads./choices.(allNames{q}).lowreward.total;
        choices.(allNames{q}).lowreward.pickedlow = choices.(allNames{q}).lowreward.pickedlow./choices.(allNames{q}).lowreward.total;
        choices.(allNames{q}).lowreward.pickedhigh = choices.(allNames{q}).lowreward.pickedhigh./choices.(allNames{q}).lowreward.total;
        choices.(allNames{q}).lowreward.correctbet = choices.(allNames{q}).lowreward.correctbet./choices.(allNames{q}).lowreward.total;
        choices.(allNames{q}).lowreward.incorrectbet = choices.(allNames{q}).lowreward.incorrectbet./choices.(allNames{q}).lowreward.total;
        choices.(allNames{q}).highreward.drewbeads = choices.(allNames{q}).highreward.drewbeads./choices.(allNames{q}).highreward.total;
        choices.(allNames{q}).highreward.pickedlow = choices.(allNames{q}).highreward.pickedlow./choices.(allNames{q}).highreward.total;
        choices.(allNames{q}).highreward.pickedhigh = choices.(allNames{q}).highreward.pickedhigh./choices.(allNames{q}).highreward.total;
        choices.(allNames{q}).highreward.correctbet = choices.(allNames{q}).highreward.correctbet./choices.(allNames{q}).highreward.total;
        choices.(allNames{q}).highreward.incorrectbet = choices.(allNames{q}).highreward.incorrectbet./choices.(allNames{q}).highreward.total;
        highdrewbeads = highdrewbeads+choices.(allNames{q}).highreward.drewbeads;
        highpickedlow = highpickedlow+choices.(allNames{q}).highreward.pickedlow;
        highpickedhigh = highpickedhigh+choices.(allNames{q}).highreward.pickedhigh;
        highcorrect = highcorrect + choices.(allNames{q}).highreward.correctbet;
        lowdrewbeads = lowdrewbeads+choices.(allNames{q}).lowreward.drewbeads;
        lowpickedlow = lowpickedlow+choices.(allNames{q}).lowreward.pickedlow;
        lowpickedhigh = lowpickedhigh+choices.(allNames{q}).lowreward.pickedhigh;
        lowcorrect = lowcorrect + choices.(allNames{q}).lowreward.correctbet;
        lowmean(q) = sum(choices.(allNames{q}).lowreward.drewbeads.*x.*choices.(allNames{q}).lowreward.total)/lownumdraws;
        highmean(q) = sum(choices.(allNames{q}).highreward.drewbeads.*x.*choices.(allNames{q}).highreward.total)/highnumdraws;
        lowsum(q) = lownumdraws/sum(choices.(allNames{q}).lowreward.total);
        highsum(q) = highnumdraws/sum(choices.(allNames{q}).highreward.total);
        if 1 % individual decisions
            subplot(6,7,q*2-1);
            title(strcat(allNames{q}, ' LowReward'),'FontSize', 20)
            hold on
            plot(x,choices.(allNames{q}).lowreward.drewbeads, 'g','LineWidth',3);
            plot(x,choices.(allNames{q}).lowreward.pickedlow, 'r','LineWidth',3);
            plot(x,choices.(allNames{q}).lowreward.pickedhigh, 'b','LineWidth',3);
            line([lowmean(q) lowmean(q)], [0 1])
            ylim([0 1]);
            subplot(6,7,q*2);
            title(strcat(allNames{q}, ' HighReward'),'FontSize', 20)
            hold on
            plot(x,choices.(allNames{q}).highreward.drewbeads, 'g','LineWidth',3);
            plot(x,choices.(allNames{q}).highreward.pickedlow, 'r','LineWidth',3);
            plot(x,choices.(allNames{q}).highreward.pickedhigh, 'b','LineWidth',3);
            line([highmean(q) highmean(q)], [0 1])
            ylim([0 1]);
        end
        if 0 % individual decisions for nodraw game
            subplot(4,4,q);
            lowmean(q) = sum(choices.(allNames{q}).lowreward.drewbeads.*x.*choices.(allNames{q}).lowreward.total)/lownumdraws;
            title(strcat(allNames{q}),'FontSize', 20)
            hold on
            plot(x,choices.(allNames{q}).lowreward.drewbeads, 'g','LineWidth',3);
            plot(x,choices.(allNames{q}).lowreward.pickedlow, 'r','LineWidth',3);
            plot(x,choices.(allNames{q}).lowreward.pickedhigh, 'b','LineWidth',3);
            line([0 0], [0 1])
            line([-2 -2], [0 1])
            line([-4 -4], [0 1])
            ylim([0 1]);
        end
    end
    if 0
        figure();
        subplot(1,2,1)
        hist(lowmean,-6:1:6)
        title('Mean LowReward','FontSize', 20)
        xlim([-6 6])
        subplot(1,2,2)
        hist(highmean,-6:1:6)
        title('Mean HighReward','FontSize', 20)
        xlim([-6 6])
    end
    highdrewbeads = highdrewbeads./length(allNames);
    highpickedlow = highpickedlow./length(allNames);
    highpickedhigh = highpickedhigh./length(allNames);
    highcorrect = highcorrect./length(allNames);
    lowdrewbeads = lowdrewbeads./length(allNames);
    lowpickedlow = lowpickedlow./length(allNames);
    lowpickedhigh = lowpickedhigh./length(allNames);
    lowcorrect = lowcorrect./length(allNames);
    if 0
        figure();
        subplot(2,1,1);
        title('Mean LowReward','FontSize', 20)
        hold on
        plot(x,lowdrewbeads, 'g','LineWidth',3);
        plot(x,lowpickedlow, 'r','LineWidth',3);
        plot(x,lowpickedhigh, 'b','LineWidth',3);
        ylim([0 1]);
        subplot(2,1,2);
        title('Mean HighReward','FontSize', 20)
        hold on
        plot(x,highdrewbeads, 'g','LineWidth',3);
        plot(x,highpickedlow, 'r','LineWidth',3);
        plot(x,highpickedhigh, 'b','LineWidth',3);
        ylim([0 1]);
    end
    
    if 0 % percentage correct
        figure();
        hold on
        plot(x,lowcorrect, 'r','LineWidth',3);
        plot(x,highcorrect, 'b','LineWidth',3);
        title('Mean correct','FontSize', 20)
        ylim([0 1]);
    end
    
    if 0 % correlation between fraction of draws(more draws = positive) and mean of the draw curve(no bias=positive).
        [r1, p1] = corr(lowmean',lowsum','type', 'Spearman','rows','complete');
        [r2, p2] = corr(highmean',highsum','type', 'Spearman','rows','complete');
    end
end

%% make regressors and write EV files
clear classes
clc
load totalData.mat
%subjects which you have scan data for
% participant NN3511 is considered a pilot for scan, and TA3514 didn't draw
% anything it seems.{'TA3514'};
list = [{'LX1803'};{'NS969'};{'QB2027'};{'OK547'};{'HU2138'}];

for m = 1: length(list)
    totalgraph = zeros(4,6700);
    for runnumber = 1:2
        for blocknumber = 1:2
            subjID = list{m};
            if runnumber == 1
                if blocknumber ==1
                    statusdata = totalData.(subjID).realScan1.block1.statusData;
                else
                    statusdata = totalData.(subjID).realScan1.block2.statusData;
                end
            else
                if blocknumber ==1
                    statusdata = totalData.(subjID).realScan2.block1.statusData;
                else
                    statusdata = totalData.(subjID).realScan2.block2.statusData;
                end
            end
            % a test to see if he/she drew at least once
            for tester = 1:100
                if ~isempty(statusdata(tester).curr_trialInfList)
                    break
                end
                if tester == 100
                    disp('this participant didnot draw anything')
                    keyboard
                end
            end
            [entirereg, fslreg] = getRegressors(statusdata);
            currentdir = cd;
            savedir = strcat('/Users/Carmenere/Documents/MATLAB/dataanalyzer for BeadsTask/beadstaskdata/',subjID, '/','EVs');
            mkdir(savedir);
            cd(savedir);
            allreg = fieldnames(fslreg);
            for q = 1:length(allreg)
                allmod = fieldnames(fslreg.(allreg{q}));
                for h = 1: length(allmod)
                    dlmwrite([subjID '_run' num2str(runnumber) '_block' num2str(blocknumber) '_' allreg{q} allmod{h} ], fslreg.(allreg{q}).(allmod{h}), ' ')
                end
            end
            cd(currentdir);
            if 0 % box car regressor calculation
                x = entirereg(:,1);
                reg1 = entirereg(:,2);
                reg2 = entirereg(:,3);
                reg3 = entirereg(:,4);
                reg4 = entirereg(:,5);
                reg5 = entirereg(:,6);
                reg6 = entirereg(:,7);
                reg7 = entirereg(:,8);
                if length(reg1)>= 6700
                    keyboard
                else
                    reg1(length(reg1)+1:6700,1) = 0;
                    reg2(length(reg2)+1:6700,1) = 0;
                    reg3(length(reg3)+1:6700,1) = 0;
                    reg4(length(reg4)+1:6700,1) = 0;
                    reg5(length(reg5)+1:6700,1) = 0;
                    reg6(length(reg6)+1:6700,1) = 0;
                    reg7(length(reg7)+1:6700,1) = 0;
                end
                temp = (1 .* reg1) + (0.9 .* reg2) + (0.8 .* reg3) + (0.7 .* reg4) + (0.6 .* reg5) + (0.5 .* reg6) + (0.4 .* reg7);
                totalgraph(2*(runnumber-1)+blocknumber,:) = totalgraph(2*(runnumber-1)+blocknumber,:) + temp';
                
                if 0 % drawing box car regressor
                    set(0,'DefaultFigureWindowStyle','docked')
                    figure();
                    subplot(7,1,1)
                    plot(x,reg1,'r')
                    xlim([0,630])
                    subplot(7,1,2)
                    plot(x,reg2,'g')
                    xlim([0,630])
                    subplot(7,1,3)
                    plot(x,reg3,'b')
                    xlim([0,630])
                    subplot(7,1,4)
                    plot(x,reg4,'r')
                    xlim([0,630])
                    subplot(7,1,5)
                    plot(x,reg5,'g')
                    xlim([0,630])
                    subplot(7,1,6)
                    plot(x,reg6,'b')
                    xlim([0,630])
                    subplot(7,1,7)
                    plot(x,reg7,'m')
                    xlim([0,630])
                end
            end
        end
    end
    if 0 %color graph for regressors
        set(0,'DefaultFigureWindowStyle','docked')
        colormap(jet)
        subplot(5,1,m), imagesc(totalgraph)
        set(gca, 'Fontsize', 16)
    end
end


