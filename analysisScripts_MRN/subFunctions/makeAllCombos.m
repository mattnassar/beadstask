function allCombosStruct=makeAllCombos(varStruct)

a=fieldnames(varStruct);

numDims=length(a);
dims=1:numDims;
data={};
% unpack all conditions and figure out how many reps you need.
for i = 1:numDims;
    eval(sprintf('data{i}=varStruct.%s;', a{i})) ;
    % make sure everything is in the first dimension.
    if size(data{i},2)>size(data{i},1);
    data{i}=data{i}';
    end
    
    numConds(i)=length(data{i});  
end

repData={};
for i = 1:numDims;
    reps=numConds;
    reps(i)=1;

    flipDims=dims;
    flipDims(i)=1;
    flipDims(1)=i;
    
    % put data in correct dimension;
    permuteData=permute(data{i}, [flipDims]);    
    matDat=repmat(permuteData, reps);
    repData{i}=matDat(:);
    eval(sprintf('allCombosStruct.%s=repData{i};', a{i})) 
end

    
    
    
    
    


