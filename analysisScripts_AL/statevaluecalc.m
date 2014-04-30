close all
clear classes
clc

limit = 30;
inc = 2;
probabilityprecision = 0;
depthlimit = 10;
calclowbound = -30;
calchighbound = 30;

UP = .60;
IC = 0.1;
HP = 70;
LP = 10;
HIP = 0;
LIP = 0;

alpha = 1; %1 is normative
b = 0; %0 is normative
cv = 0; %0 is normative

x=1;
tic
func1 = nan(1,2*limit/inc+1);
func2 = nan(1,2*limit/inc+1);
func3 = nan(1,2*limit/inc+1);
for BD = -limit:inc:limit
    [func1(x), func2(x), func3(x)] = calculator(BD, UP, IC, HP, LP, HIP, LIP, limit, 1, alpha, b, cv, probabilityprecision, depthlimit, 0, calclowbound, calchighbound);
    x/length(-limit:inc:limit)
    x=x+1;
end
toc
hold on
BD = -limit:inc:limit;
maxfunc = max(func1, func2);
area(BD, 50000*(func2> func3 & func2>func1), 'FaceColor', [1 0.7 0.7], 'LineWidth', 4, 'EdgeColor', [1 1 1])
area(BD, 50000*(func1> func3 & func1>func2), 'FaceColor', [0.7 0.7 1], 'LineWidth', 4, 'EdgeColor', [1 1 1])
area(BD, 50000*(func3> maxfunc), 'FaceColor', [0.7 1 0.7], 'LineWidth', 4, 'EdgeColor', [1 1 1])
difffunc = func3 -maxfunc;

difffunc = func3 - maxfunc;
title('Optimal Decision Model', 'FontSize', 30)
xlabel('Bead Difference', 'FontSize', 24)
ylabel('Value', 'FontSize', 24)

line([-limit limit]',[0 0]','LineWidth', 2, 'LineStyle', '-.', 'Color', [0 0 0])
x = (find(difffunc == max(difffunc))-1) *inc - limit;
if length(x) == 1 && isfinite(x)
    line([x x]',[0 70]','LineWidth', 2, 'LineStyle', ':', 'Color', [0.4 0.4 0.4])
end
line([0 0]',[-70 70]','LineWidth', 2, 'LineStyle', '-.', 'Color', [0 0 0])

difffunc(find(difffunc < 0)) = 0;
infVal = unitNorm(difffunc);
%plot(BD,70.*infVal, 'm','LineWidth', 2)
plot(BD, func3, 'g', 'LineWidth', 8)
plot(BD, func1, 'b', 'LineWidth', 8)
plot(BD, func2, 'r', 'LineWidth', 8)
xlim([-15 15])
ylim([0 70])