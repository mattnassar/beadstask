function [allData] = indivFullDataLoader(dirandfile)
everything=dir(fullfile(dirandfile));
fn={everything.name};
allData = struct;
cd(dirandfile);
flag = 0;
for q = 1:length(fn)
    if length(fn{q}) >20 && ~strcmp(fn{q}(end-14:end),'topsDataLog.mat')
        if ~isempty(strfind(fn{q},'noBead'))
            if flag == 0
                allData.noBead.block1 = load(fn{q});
                flag = 1;
            else
                allData.noBead.block2 = load(fn{q});
            end
        elseif ~isempty(strfind(fn{q},'symm'))
            allData.symm.block1 = load(fn{q});
        else
            blocknumber = fn{q}(sum([diff(strfind(fn{q},'_'))==2 0].*strfind(fn{q},'_'))+1);
            % this looks for a place in the string where two underbars have
            % one character in between. That character is the blocknumber
            if ~isempty(strfind(fn{q},'noDraw'))
                eval(sprintf('allData.noDraw.block%s = load(fn{q});',blocknumber));
            elseif ~isempty(strfind(fn{q},'game1'))
                eval(sprintf('allData.game1.block%s = load(fn{q});',blocknumber));
            elseif ~isempty(strfind(fn{q},'scan'))
                eval(sprintf('allData.scan.block%s = load(fn{q});',blocknumber));
            else
                realScanversion = fn{q}(strfind(fn{q},'realScannerGame')+15);
                eval(sprintf('allData.realScan%s.block%s = load(fn{q});',realScanversion,blocknumber));
            end
        end
    end
end
[allData] = orderfields(allData, {'noBead','noDraw', 'symm', 'game1', 'scan','realScan0','realScan1','realScan2'});



