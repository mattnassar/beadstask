clear classes
clc
cd /Users/Carmenere/Documents/MATLAB
load choices.mat % created by shell.m which analyzed drawing behavior
pilots = [{'BJ2376'};{'BX2031'};{'DT2005'};{'NC3509'};{'SQ3510'}]; %%these five were pilots
for index1 = 1:length(pilots) %removing pilots from analysis
    choices = rmfield(choices, pilots{index1});
end
allNames = fieldnames(choices);
x=-12:2:12;
%% draw curve
figure();
for index1 = 1:length(allNames)
    lownumdraws = sum(choices.(allNames{index1}).lowreward.drewbeads);
    BD = -10:2:10;
    centerofcurve(index1) = sum(choices.(allNames{index1}).lowreward.drewbeads.*BD)/lownumdraws;
    choices.(allNames{index1}) = choices.(allNames{index1}).lowreward.drewbeads./choices.(allNames{index1}).lowreward.total;
    y = [0 choices.(allNames{index1}) 0];
    subplot(4,4,index1);
    hold on
    plot(x,y, 'g','LineWidth',3);
    line([0 0]',[0 1]','Color',[0 0 0], 'LineWidth', 2, 'LineStyle', '-.')
    plot(x,y, 'g','LineWidth',3);
    line([centerofcurve(index1) centerofcurve(index1)], [0 1])
    title(strcat(allNames{index1}),'FontSize',10);
    xlim([-12 12])
    ylim([0 1])
    hold off
end


%% power utility fitting
figure();
beginningalpha = nan(1,length(allNames));
for index1 = 1:length(allNames)
    beginningalpha(index1) = centerofcurve(index1)*(log(3/2)/log(1/7));
    if beginningalpha(index1) < 0
        beginningalpha(index1) = 0.1;
    end
    if beginningalpha(index1) > 2
        beginningalpha(index1) = 2;
    end
    beginningalpha(index1)
    y = [0 choices.(allNames{index1}) 0];
    subplot(4,4,index1);
    hold on
    plot(x,y, 'g','LineWidth',3);
    line([0 0]',[0 1]','Color',[0 0 0], 'LineWidth', 2, 'LineStyle', '-.')
    plot(x,y, 'g','LineWidth',3);
    
    params = [beginningalpha(index1),1];
    model = @(params,x) fitcurve(params(1),0 , 0, params(2));  %fitcurve(alpha, b, cv, invTemp, gain)
    [fitparams, resnorm] = lsqcurvefit(model,params,x, y,[0 0],[2 Inf]);
    plot(x, fitcurve(fitparams(1),0, 0, fitparams(2)),'b');
    alpha = strcat('a = ',num2str(roundn(fitparams(1),-2)));
    inverse = strcat('iT = ', num2str(round(fitparams(2))));
    res = strcat('res = ', num2str(roundn(resnorm, -2)));
    title(strcat(allNames{index1},',',alpha,',',inverse,',',res),'FontSize',10);
    xlim([-12 12])
    ylim([0 1])
    hold off
end

%% normative risk-return
figure();
beginningb = nan(1,length(allNames));
for index1 = 1:length(allNames)
    beginningb(index1) = (7*(3/2)^centerofcurve(index1) + 6 - (2/3)^centerofcurve(index1))/480;
    beginningb(index1)
    y = [0 choices.(allNames{index1}) 0];
    subplot(4,4,index1);
    hold on
    plot(x,y, 'g','LineWidth',3);
    line([0 0]',[0 1]','Color',[0 0 0], 'LineWidth', 2, 'LineStyle', '-.')
    plot(x,y, 'g','LineWidth',3);
    
    params = [beginningb(index1),1];
    model = @(params,x) fitcurve(1,params(1), 0, params(2));  %fitcurve(alpha, b, cv, invTemp)
    [fitparams, resnorm] = lsqcurvefit(model,params,x, y,[0 0],[Inf Inf]);
    plot(x, fitcurve(1,fitparams(1), 0, fitparams(2)),'b');
    alpha = strcat('b = ',num2str(roundn(fitparams(1),-4)));
    inverse = strcat('iT = ', num2str(round(fitparams(2))));
    res = strcat('res = ', num2str(roundn(resnorm, -2)));
    title(strcat(allNames{index1},',',alpha,',',inverse,',',res),'FontSize',10);
    xlim([-12 12])
    ylim([0 1])
    hold off
end

%% Weber's risk-return
figure();
beginningcv = nan(1,length(allNames));
for index1 = 1:length(allNames)
    z = centerofcurve(index1);
    beginningcv(index1) = (6^(z/2)*(70*3^z-10*2^z))/(2^(2*z)-3^(2*z));
    beginningcv(index1)
    y = [0 choices.(allNames{index1}) 0];
    subplot(4,4,index1);
    hold on
    plot(x,y, 'g','LineWidth',3);
    line([0 0]',[0 1]','Color',[0 0 0], 'LineWidth', 2, 'LineStyle', '-.')
    plot(x,y, 'g','LineWidth',3);
    
    params = [1,1];
    model = @(params,x) fitcurve(1,0,params(1), params(2));  %fitcurve(alpha, b, cv, invTemp)
    [fitparams, resnorm] = lsqcurvefit(model,params,x, y,[0 0],[Inf Inf]);
    plot(x, fitcurve(1,0,fitparams(1), fitparams(2)),'b');
    alpha = strcat('Wb = ',num2str(roundn(fitparams(1),-2)));
    inverse = strcat('iT = ', num2str(round(fitparams(2))));
    res = strcat('res = ', num2str(roundn(resnorm, -2)));
    title(strcat(allNames{index1},',',alpha,',',inverse,',',res),'FontSize',10);
    xlim([-12 12])
    ylim([0 1])
    hold off
end

%% Symmetric game
figure();
for index1 = 1:length(allNames)
    index1
    y = [0 choices.(allNames{index1}) 0];
    subplot(4,4,index1);
    hold on
    plot(x,y, 'g','LineWidth',3);
    line([0 0]',[0 1]','Color',[0 0 0], 'LineWidth', 2, 'LineStyle', '-.')
    plot(x,y, 'g','LineWidth',3);
    
    params = [70,1];
    model = @(params,x) fitcurve(1,0,0, params(2),params(1),params(1));  %fitcurve(alpha, b, cv, invTemp, HP, LP)
    [fitparams, resnorm] = lsqcurvefit(model,params,x, y,[0 0],[70 Inf]);
    plot(x, fitcurve(1,0,0,fitparams(2), fitparams(1), fitparams(1)),'b');
    alpha = strcat('Pay = ',num2str(round(fitparams(1))));
    inverse = strcat('iT = ', num2str(round(fitparams(2))));
    res = strcat('res = ', num2str(roundn(resnorm, -2)));
    title(strcat(allNames{index1},',',alpha,',',inverse,',',res),'FontSize',10);
    xlim([-12 12])
    ylim([0 1])
    hold off
end
