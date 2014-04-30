close all
clear classes
clc
cd '/Users/Carmenere/Documents/MATLAB/dataanalyzer for BeadsTask'
participantlist = [{'BZ1909'};{'FI3535'};{'WE3522'};{'IL3521'};{'NG977'};{'TX1901'};{'YQ3524'};{'ML494'}];

for m = 1:length(participantlist)
    subjID = participantlist{m,:};
    [allData] = revisedbeadTaskDataLoader(subjID);
    eval(sprintf('totalData.%s=allData', subjID));
end
cd '/Users/Carmenere/Documents/MATLAB/dataanalyzer for BeadsTask'
save totalDataforrevisedtask.mat totalData

%% individual task graph
clear classes
clc
load totalDataforrevisedtask.mat
allNames = fieldnames(totalData);
shell.bethigh = zeros(1,11);
shell.betlow = zeros(1,11);
shell.draw = zeros(1,11);
shell.total = zeros(1,11);
for m = 1:length(allNames) % subject level
    choices.nodraw60 = shell;
    choices.nodraw70 = shell;
    choices.nodraw75 = shell;
    choices.game1low = shell;
    choices.game1high = shell;
    choices.game2low = shell;
    choices.game2high = shell;
    choices.scanlow = shell;
    choices.scanhigh = shell;
    fn = fieldnames(totalData.(allNames{m}));
    for n = 1:length(fn)
        statusData = [totalData.(allNames{m}).(fn{n}).block1.statusData;totalData.(allNames{m}).(fn{n}).block2.statusData];
        for q = 1:length(statusData)
            if ~isempty(statusData(q).curr_rewCorrLeft) % if this is not an empty shell
                if statusData(q).curr_rewCorrLeft > statusData(q).curr_rewCorrRight
                    highside = 1;
                else
                    highside = 2;
                end
                if highside == 1
                    BD = (statusData(q).curr_start_tokensLeft - statusData(q).curr_start_tokensRight)/2+6;
                else
                    BD = (statusData(q).curr_start_tokensRight - statusData(q).curr_start_tokensLeft)/2+6;
                end
                firstchoice = 0; %1 is bet low, 2 is draw, 3 is bet high, 0 is timed out
                if isempty(statusData(q).curr_trialInfList) % didn't draw
                    if statusData(q).isGoodTrial % choice was made
                        if statusData(q).curr_choice == highside % bet high
                            firstchoice = 3;
                        else % bet low
                            firstchoice = 1;
                        end
                    end
                else % drew beads
                    firstchoice = 2;
                end
                string = '';
                switch max(statusData(q).curr_rewCorrLeft, statusData(q).curr_rewCorrRight)
                    case 60 % nodraw60
                        string = 'nodraw60';
                    case 75 % nodraw75
                        string = 'nodraw75';
                    case 105 % game 1 high pay
                        string = 'game1high';
                    case 7 % game 2 low pay
                        string = 'game2low';
                    case 170 % scan high pay
                        string = 'scanhigh';
                    case 70 % nodraw70, game1 low pay, game2 high pay, scan low pay
                        switch n
                            case 2 % nodraw70
                                string = 'nodraw70';
                            case 4 % game1 low pay
                                string = 'game1low';
                            case 5 % game2 high pay
                                string = 'game2high';
                            case 6 % scan low pay
                                string = 'scanlow';
                        end
                end
                switch firstchoice
                    case 0 % timed out
                        [m,n,q]
                    case 1 % bet low
                        choices.(string).betlow(BD) = choices.(string).betlow(BD)+1;
                        choices.(string).total(BD) = choices.(string).total(BD)+1;
                    case 2 % draw
                        choices.(string).draw(BD) = choices.(string).draw(BD)+1;
                        choices.(string).total(BD) = choices.(string).total(BD)+1;
                    case 3 % bet high
                        choices.(string).bethigh(BD) = choices.(string).bethigh(BD) +1;
                        choices.(string).total(BD) = choices.(string).total(BD)+1;
                    otherwise
                        keyboard
                end
            end
        end
    end
    figure();
    subplot(9,1,1)
    hold on
    plot(-10:2:10, choices.nodraw60.betlow, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.nodraw60.bethigh, 'b', 'LineWidth', 3)
    line([-2.7 -2.7], [0 8])
    ylim([0 8])
    title(strcat('Subject : ', allNames{m}, ' _ nodraw _ 60:20 _ optimal = -2.7'), 'FontSize', 16)
    
    subplot(9,1,2)
    hold on
    plot(-10:2:10, choices.nodraw70.betlow, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.nodraw70.bethigh, 'b', 'LineWidth', 3)
    line([-4.8 -4.8], [0 8])
    ylim([0 8])
    title(strcat('Subject : ', allNames{m}, ' _ nodraw _ 70:10 _ optimal = -4.8'), 'FontSize', 16)
    
    subplot(9,1,3)
    hold on
    plot(-10:2:10, choices.nodraw75.betlow, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.nodraw75.bethigh, 'b', 'LineWidth', 3)
    line([-6.7 -6.7], [0 8])
    ylim([0 8])
    title(strcat('Subject : ', allNames{m}, ' _ nodraw _ 75:5 _ optimal = -6.7'), 'FontSize', 16)
    
    subplot(9,1,4)
    hold on
    plot(-10:2:10, choices.game1low.betlow, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.game1low.bethigh, 'b', 'LineWidth', 3)
    plot(-10:2:10, choices.game1low.draw, 'g', 'LineWidth', 3)
    line([-4.8 -4.8], [0 8])
    ylim([0 4])
    title(strcat('Subject : ', allNames{m}, ' _ game1 _ 70:10 _ optimal = -4.8'), 'FontSize', 16)
    
    subplot(9,1,5)
    hold on
    plot(-10:2:10, choices.game1high.betlow, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.game1high.bethigh, 'b', 'LineWidth', 3)
    plot(-10:2:10, choices.game1high.draw, 'g', 'LineWidth', 3)
    line([-4.8 -4.8], [0 8])
    ylim([0 4])
    title(strcat('Subject : ', allNames{m}, ' _ game1 _ 105:45 _ optimal = -4.8'), 'FontSize', 16)
    
    subplot(9,1,6)
    hold on
    plot(-10:2:10, choices.game2low.betlow, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.game2low.bethigh, 'b', 'LineWidth', 3)
    plot(-10:2:10, choices.game2low.draw, 'g', 'LineWidth', 3)
    line([-4.8 -4.8], [0 8])
    ylim([0 4])
    title(strcat('Subject : ', allNames{m}, ' _ game2 _ 7:1 _ optimal = -4.8'), 'FontSize', 16)
    
    subplot(9,1,7)
    hold on
    plot(-10:2:10, choices.game2high.betlow, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.game2high.bethigh, 'b', 'LineWidth', 3)
    plot(-10:2:10, choices.game2high.draw, 'g', 'LineWidth', 3)
    line([-4.8 -4.8], [0 8])
    ylim([0 4])
    title(strcat('Subject : ', allNames{m}, ' _ game2 _ 70:10 _ optimal = -4.8'), 'FontSize', 16)
    
    subplot(9,1,8)
    hold on
    plot(-10:2:10, choices.scanlow.betlow./choices.scanlow.total, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.scanlow.bethigh./choices.scanlow.total, 'b', 'LineWidth', 3)
    plot(-10:2:10, choices.scanlow.draw./choices.scanlow.total, 'g', 'LineWidth', 3)
    line([-4.8 -4.8], [0 8])
    ylim([0 1])
    title(strcat('Subject : ', allNames{m}, ' _ scan _ 70:10 _ optimal = -4.8'), 'FontSize', 16)
    
    subplot(9,1,9)
    hold on
    plot(-10:2:10, choices.scanhigh.betlow./choices.scanhigh.total, 'r', 'LineWidth', 3)
    plot(-10:2:10, choices.scanhigh.bethigh./choices.scanhigh.total, 'b', 'LineWidth', 3)
    plot(-10:2:10, choices.scanhigh.draw./choices.scanhigh.total, 'g', 'LineWidth', 3)
    line([-4.8 -4.8], [0 8])
    ylim([0 1])
    title(strcat('Subject : ', allNames{m}, ' _ scan _ 170:110 _ optimal = -4.8'), 'FontSize', 16)
end
