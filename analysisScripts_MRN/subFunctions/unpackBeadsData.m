function [dataStruct]=unpackBeadsData(statusData)
% keyboard
% find the first good trial.
k=1;
while (isempty(statusData(k).isGoodTrial)| ~statusData(k).isGoodTrial )&k<100
    k=k+1;
end

if k<99
    
    ll=length([statusData.curr_maxUrnProb]);
    
    allNames=fieldnames(statusData);
    dataStruct.curr_trialInfList={statusData.curr_trialInfList};
    for j = 1:length(allNames)
        % get data from a field
        if (isnumeric(eval(sprintf('statusData(%d).%s', k, char(allNames(j)))))||...
                islogical(eval(sprintf('statusData(%d).%s', k, char(allNames(j))))))...
                && (length((eval(sprintf('statusData(%d).%s', k, char(allNames(j))))))==1)...
                && isempty(strmatch(char(allNames(j)), {'extraDrawTime', 'curr_trialInfList'})) ;
            eval(sprintf('dataStruct.%s=cat(1, statusData.%s);', char(allNames(j)), char(allNames(j))));
            
            %% if the length is wrong, fix it with nans.
            lDiff=ll - length(eval(sprintf('dataStruct.%s', char(allNames(j)))));
            if lDiff>0
                disp('problem with lengths, padding with nans')
                eval(sprintf('dataStruct.%s=cat(1, nan(lDiff, 1), dataStruct.%s);', char(allNames(j)), char(allNames(j))))
            end
        elseif (isnumeric(eval(sprintf('statusData(%d).%s', k, char(allNames(j)))))||...
                islogical(eval(sprintf('statusData(%d).%s', k, char(allNames(j))))))...
                && length((eval(sprintf('statusData(%d).%s', k, char(allNames(j))))))>1;
            eval(sprintf('dataStruct.%s={statusData.%s};', char(allNames(j)), char(allNames(j))));
        end
    end
    dataStruct=straightStruct(dataStruct);
    
else
    dataStruct=[]
    input('missing subject data')
end
