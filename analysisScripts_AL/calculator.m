function [bethigh, betlow, draw] = calculator(BD, UP, IC, HP, LP, HIP, LIP, limit, currprob, alpha, b, cv, probabilityprecision, depthlimit, depth, calclowbound, calchighbound)
depth = depth +1;
probHU = 1/(1+((1-UP)/UP)^BD);
probLU = 1/(1+((1-UP)/UP)^(-BD));
probHB = probHU*UP + probLU * (1-UP);
probLB = probLU*UP + probHU * (1-UP);
EXPbethigh = (probHU*HP^alpha + probLU*HIP^alpha);
EXPbetlow = (probLU*LP^alpha + probHU*LIP^alpha);
if depth > 1
    bethigh = EXPbethigh; %return EV
    betlow = EXPbetlow; %return EV
else
    if depth == 1
        bethighvar = probHU *(HP-EXPbethigh)^2 + probLU*(HIP-EXPbethigh)^2;
        bethigh = EXPbethigh - b*bethighvar - cv *sqrt(bethighvar)/EXPbethigh; %return WTP
        betlowvar = probLU *(LP-EXPbetlow)^2 + probHU*(LIP-EXPbetlow)^2;
        betlow = EXPbetlow - b*betlowvar - cv * sqrt(betlowvar)/EXPbetlow; %return WTP
    else
        keyboard
    end
end
currprobHB = currprob * probHB;
currprobLB = currprob * probLB;
if BD == -(limit)
    draw = max(LP,HIP)-IC;
elseif BD == limit
    draw = max(HP,LIP)-IC;
else
    if currprobHB > probabilityprecision && depth <= depthlimit
        if BD+1 >= calchighbound
            PD = BD+1;
            newprobHU = 1/(1+((1-UP)/UP)^PD);
            newprobLU = 1/(1+((1-UP)/UP)^(-PD));
            newEXPbethigh = (newprobHU*HP^alpha + newprobLU*HIP^alpha);
            newEXPbethighvar = newprobHU *(HP-newEXPbethigh)^2 + newprobLU*(HIP-newEXPbethigh)^2;
            SBDplus1 = newEXPbethigh - b*(newEXPbethighvar) - cv *sqrt(newEXPbethighvar)/newEXPbethigh;
        else
            [output1, output2, output3] = calculator(BD+1, UP, IC, HP, LP, HIP, LIP, limit,currprobHB, alpha, b, cv, probabilityprecision, depthlimit, depth, calclowbound, calchighbound);
            SBDplus1 = max(max(output1, output2), output3);
        end
    else
        SBDplus1 = 0;
    end
    if currprobLB > probabilityprecision && depth <= depthlimit
        if BD-1 <= calclowbound
            PD = BD-1;
            newprobHU = 1/(1+((1-UP)/UP)^PD);
            newprobLU = 1/(1+((1-UP)/UP)^(-PD));
            newEXPbetlow = (newprobLU*LP^alpha + newprobHU*LIP^alpha);
            newEXPbetlowvar = newprobLU *(LP-newEXPbetlow)^2 + newprobHU*(LIP-newEXPbetlow)^2;
            SBDneg1 = newEXPbetlow - b*(newEXPbetlowvar) - cv * sqrt(newEXPbetlowvar)/newEXPbetlow;
        else
            [output1, output2, output3] = calculator(BD-1, UP, IC, HP, LP, HIP, LIP, limit,currprobLB, alpha, b, cv, probabilityprecision, depthlimit, depth, calclowbound, calchighbound);
            SBDneg1 = max(max(output1, output2), output3);
        end
    else
        SBDneg1 = 0;
    end
    
    EXPdraw = probHB*SBDplus1+probLB*SBDneg1-IC^alpha;
    if depth > 1
        draw = EXPdraw; %return EV
    else
        if depth == 1
            drawvar = probHB * (SBDplus1-EXPdraw)^2 + probLB * (SBDneg1 - EXPdraw)^2;
            draw = EXPdraw - b*drawvar - cv*sqrt(drawvar)/EXPdraw; %return WTP
        else
            keyboard
        end
    end
end


