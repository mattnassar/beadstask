function [draw] = todrawornot(beaddiff)
if beaddiff >= 8 | beaddiff <= -12
    draw = 0;
else
%     if this code line is reached, beaddiff is an integer in the range of
%     -11<= beaddiff<= 7
    interpX=  -11:1:7;
    F=[0.0893    0.1786    0.2445    0.3103    0.3649    0.4194 ...   
   0.4855    0.5517    0.4828    0.4138    0.2931    0.1724    0.1668 ...
   0.1613    0.1306    0.1000    0.0969    0.0938    0.0469]; 
    
    if rand(1) <= F(interpX==beaddiff)
        draw = 1;
    else
        draw = 0;
    end
end

