function [allData] = revisedbeadTaskDataLoader(subjID)
dirandfile = strcat('/Users/Carmenere/Documents/MATLAB/dataanalyzer for BeadsTask/revisedbeadsdata/',subjID);
everything=dir(fullfile(dirandfile));
fn={everything.name};
allData = struct;
cd(dirandfile);
for q = 1:length(fn)
    if length(fn{q}) >20
        tempstring = fn{q};
        matcher = tempstring(length(subjID)+2:length(subjID)+6);
        if strcmp(matcher, 'game1')
            matcher = tempstring(length(subjID)+8);
            if strcmp(matcher, '1')
                % game 1 block 1
                allData.game1.block1 = load(tempstring);
            elseif strcmp(matcher, '2')
                % game 1 block 2
                allData.game1.block2 = load(tempstring);
            else
                keyboard
            end
        elseif strcmp(matcher, 'game2')
            matcher = tempstring(length(subjID)+8);
            if strcmp(matcher, '1')
                % game 2 block 1
                allData.game2.block1 = load(tempstring);
            elseif strcmp(matcher, '2')
                % game 2 block 2
                allData.game2.block2 = load(tempstring);
            else
                keyboard
            end
        elseif strcmp(matcher, 'noDra')
            matcher = tempstring(length(subjID)+26:length(subjID)+27);
            if strcmp(matcher, '75')
                matcher = tempstring(length(subjID)+19);
                if strcmp(matcher, '1')
                    % no draw block 1
                    allData.nodraw75.block1 = load(tempstring);
                elseif strcmp(matcher, '2')
                    % no draw block 2
                    allData.nodraw75.block2 = load(tempstring);
                else
                    keyboard
                end
            elseif strcmp(matcher, '70')
                matcher = tempstring(length(subjID)+19);
                if strcmp(matcher, '1')
                    % no draw block 1
                    allData.nodraw70.block1 = load(tempstring);
                elseif strcmp(matcher, '2')
                    % no draw block 2
                    allData.nodraw70.block2 = load(tempstring);
                else
                    keyboard
                end
            elseif strcmp(matcher, '60')
                matcher = tempstring(length(subjID)+19);
                if strcmp(matcher, '1')
                    % no draw block 1
                    allData.nodraw60.block1 = load(tempstring);
                elseif strcmp(matcher, '2')
                    % no draw block 2
                    allData.nodraw60.block2 = load(tempstring);
                else
                    keyboard
                end
            else
                keyboard
            end
        elseif strcmp(matcher, 'scann')
            matcher = tempstring(length(subjID)+16);
            if strcmp(matcher, '1')
                % scan block 1
                allData.scan.block1 = load(tempstring);
            elseif strcmp(matcher, '2')
                % scan block 2
                allData.scan.block2 = load(tempstring);
            else
                keyboard
            end
        else
            keyboard
        end
    end
end

[allData] = orderfields(allData, {'nodraw60', 'nodraw70', 'nodraw75', 'game1', 'game2', 'scan'});