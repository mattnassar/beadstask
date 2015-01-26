function [allData, p, t] = beadTaskDataLoader_mrn(dirandfile)

[~, subjID]=fileparts(dirandfile)


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
        elseif strcmp(matcher, 'noDraw')
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
        elseif strcmp(matcher, 'realS')
            label='realScannerGame';
            labelStart=strfind(tempstring, label);
            labelLength=length(label);
            eval(sprintf('allData.%s.block%s = load(tempstring)', tempstring(labelStart:labelStart+labelLength), tempstring(labelStart+labelLength+3)));                 
        elseif strcmp(matcher, 'symm_')
            % symm condition
            allData.symm.block1 = load(tempstring);
        else
            keyboard
        end
    end
end
%[allData] = orderfields(allData, {'nodraw', 'symm', 'game1', 'game2', 'scan'});

if nargout>1
%% analyze the fraction of bead draws
x = -10:2:10;
p = zeros(1,21);
t = zeros(1,21);
for m = 1:2
    if m == 1
        statusData = allData.scan.block1.statusData;
    else
        statusData = allData.scan.block2.statusData;
    end
    for q= 1:length(statusData)
        if statusData(q).curr_rewCorrLeft > statusData(q).curr_rewCorrRight
            % left option has higher payout
            temp = statusData(q).curr_start_tokensLeft - statusData(q).curr_start_tokensRight;
            % temp is the number of beads favoring the high payout option
            t(temp+11) = t(temp+11)+1;
            if ~isempty(statusData(q).curr_trialInfList)
                % the participant drew some beads
                p(temp+11) = p(temp+11)+1;
            end
        elseif statusData(q).curr_rewCorrLeft < statusData(q).curr_rewCorrRight
            % right option has higer payout
            temp = statusData(q).curr_start_tokensRight - statusData(q).curr_start_tokensLeft;
            % temp is the number of beads favoring the high payout option
            t(temp+11) = t(temp+11)+1;
            if ~isempty(statusData(q).curr_trialInfList)
                % the participant drew some beads
                p(temp+11) = p(temp+11)+1;
            end
        elseif statusData(q).curr_rewCorrLeft == statusData(q).curr_rewCorrRight
            % equivalent payout
            keyboard
        elseif isempty(statusData(q).curr_rewCorrLeft) || isempty(statusData(q).curr_rewCorrRight)
            % empty cells
        else
            %unexpected error
            keyboard
        end
    end
end
p = p(1:2:21);
t = t(1:2:21);

end



