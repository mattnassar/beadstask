function [allData] = beadTaskDataLoader(subjID)
dirandfile = strcat('/Users/Carmenere/Documents/MATLAB/dataanalyzer for BeadsTask/beadstaskdata/',subjID);
everything=dir(fullfile(dirandfile));
fn={everything.name};
allData = struct;
cd(dirandfile);
flagforscan = 0;
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
            matcher = tempstring(length(subjID)+19);
            if strcmp(matcher, '1')
                % no draw block 1
                allData.nodraw.block1 = load(tempstring);
            elseif strcmp(matcher, '2')
                % no draw block 2
                allData.nodraw.block2 = load(tempstring);
            else
                keyboard
            end
        elseif strcmp(matcher, 'realS')
            flagforscan = 1;
            matcher = tempstring(length(subjID)+17);
            if strcmp(matcher, '1')
                % realscan game1
                matcher = tempstring(length(subjID)+20);
                if strcmp(matcher, '1')
                    allData.realScan1.block1 = load(tempstring);
                elseif strcmp(matcher, '2')
                    allData.realScan1.block2 = load(tempstring);
                else
                    keyboard
                end
            elseif strcmp(matcher, '2')
                % realscan game2
                matcher = tempstring(length(subjID)+20);
                if strcmp(matcher, '1')
                    allData.realScan2.block1 = load(tempstring);
                elseif strcmp(matcher, '2')
                    allData.realScan2.block2 = load(tempstring);
                else
                    keyboard
                end
            elseif strcmp(matcher, '0')
                % realscan game0
                matcher = tempstring(length(subjID)+20);
                if strcmp(matcher, '1')
                    allData.realScan0.block1 = load(tempstring);
                elseif strcmp(matcher, '2')
                    allData.realScan0.block2 = load(tempstring);
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
        elseif strcmp(matcher, 'symm_')
            % symm condition
            allData.symm.block1 = load(tempstring);
        else
            keyboard
        end
    end
end
if flagforscan == 1
    [allData] = orderfields(allData, {'nodraw', 'symm', 'game1', 'game2', 'scan', 'realScan0', 'realScan1', 'realScan2'});
else
    [allData] = orderfields(allData, {'nodraw', 'symm', 'game1', 'game2', 'scan'});
end