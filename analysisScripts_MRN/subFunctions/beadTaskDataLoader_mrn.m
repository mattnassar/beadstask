function [allData, p, t] = beadTaskDataLoader_mrn(dirandfile)

[~, subjID]=fileparts(dirandfile);


everything=dir(fullfile(dirandfile));
%keyboard
fn={everything.name};
allData = struct;
cd(dirandfile);


fileList={'noDraw_noBead_',   'noDraw_asym_cond_', 'symm_cond', 'game1_1', ...
    'game1_2','scannerGame1__1', 'scannerGame1__2',  ...
    'realScannerGame0__1', 'realScannerGame0__2', ...
    'realScannerGame1', 'realScannerGame2'};


for q = 1:length(fn)
    if length(fn{q}) >20
        
        if ~isempty(strfind(fn{q}, fileList{1}))
            if exist('allData')&&isfield(allData, 'noDraw_noBead')&&isfield(allData.noDraw_noBead, 'block1')
                allData.noDraw_noBead.block2 = load(fn{q}); % if there is already a block 1, then this is block 2. 
            else
                allData.noDraw_noBead.block1 = load(fn{q}); % first block first...
            end
        elseif ~isempty(strfind(fn{q}, fileList{2}))
            disp('no draw!')
            if exist('allData')&&   isfield(allData, 'noDraw_asym')  && isfield(allData.noDraw_asym, 'block1')
                allData.noDraw_asym.block2 = load(fn{q}); % if there is already a block 1, then this is block 2. 
            else
                allData.noDraw_asym.block1 = load(fn{q}); % first block first...
            end
        elseif ~isempty(strfind(fn{q}, fileList{3}))
                allData.game1.symm_cond = load(fn{q});
        elseif ~isempty(strfind(fn{q}, fileList{4}))
                allData.game1.block1 = load(fn{q});
        elseif ~isempty(strfind(fn{q}, fileList{5}))
                allData.game1.block2 = load(fn{q});
        elseif ~isempty(strfind(fn{q}, fileList{6}))
                allData.scanWarmup.block1=load(fn{q});
        elseif ~isempty(strfind(fn{q}, fileList{7}))
                allData.scanWarmup.block2=load(fn{q});
        elseif ~isempty(strfind(fn{q}, fileList{8}))
                allData.scanerGame.block1=load(fn{q});
        elseif ~isempty(strfind(fn{q}, fileList{9}))
                allData.scanerGame.block2=load(fn{q});       
        elseif ~isempty(strfind(fn{q}, fileList{10})) | ~isempty(strfind(fn{q}, fileList{11}))
            label='realScannerGame';
            labelStart=strfind(fn{q}, label);
            labelLength=length(label);
            eval(sprintf('allData.%s.block%s = load(fn{q})', fn{q}(labelStart:labelStart+labelLength), fn{q}(labelStart+labelLength+3)));                 
        else
            %keyboard
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



