%% Arthur's script
clear classes
clc

cd '/Users/Carmenere/Documents/ANALYSIS/Notes, Data and Scripts/SCRIPTS/BeadsTask/analysisScripts_AL/designfiles'

subjectIDs = [{'CC0528'},{'FL1136'},{'FT3594'},{'IL3520'},{'LL3555'},{'NS0950'},{'OF3592'},{'TK3556'},{'TQ3543'},{'TQ3600'}];
runnums = [{'r01'},{'r02'},{'r03'},{'r04'}];
numberofEVs = 9;
numberofRuns = 4;
oldRunText='r01';
oldSubjectText='TQ3543';
evFormat=[oldSubjectText '_realScannerGame1__1_'];

fileID = fopen('general_design.fsf');
file=fread(fileID, 'uint8=>char')';
fclose(fileID);
EVindex = regexp(file, evFormat);
subjectindex = regexp(file,oldSubjectText);


%change ev filenames, outputdirectory
for i = 1:length(subjectIDs)
    person = subjectIDs{i};
    for k = 1:numberofRuns
        if k <3
            gamenum = '1';
        else
            gamenum = '2';
        end
        if mod(k,2) == 0
            blocknum = '2';
        else
            blocknum = '1';
        end
        %changing EV filenames
        for j = 1:numberofEVs
            file(EVindex(j):EVindex(j)+5) = person;
            file(EVindex(j)+22) = gamenum;
            file(EVindex(j)+25) = blocknum;
        end
        
        %changing output directory
        file(subjectindex(1)+7:subjectindex(1)+9) = runnums{k};
        
        %changing functional files
        file(subjectindex(2)+7:subjectindex(2)+9) = runnums{k};
        
        for q = 1:length(subjectindex)
            file(subjectindex(q):subjectindex(q)+5) = person;
        end
        
        newDF=[person runnums{k} '_design.fsf'];
        newFID=fopen(newDF, 'w+');
        COUNT=fwrite(newFID, file);
    end
end

%change output directory

%% BIG SCRIPT TO AUTOMATE fMRI ANALYSIS!!!!

% a work in progress by MRN, started 1/07/13.  Goal is to make this script
% a template to use for future fMRI studies. The more achievable goal is to
% do a basic GLM analysis of the change-point dataset collected by JFM and
% MRN.

% the script relies on FSL tools and is meant to be run from the cluster,
% where it will operate FSL through the unix command line.

% Really run pre-processing?
reallyRun=1;
%% Setup

% specify directories where data are located


baseDir='/data/jet/kable/Arthur/';
rawDatDir='/data/jet/kable/Arthur/beadsFuncData/';
anatDir='/data/jet/kable/Arthur/beadsAnatomical/';
outputDir=[baseDir 'beads_lowLev_2-16-15_MRN'];
templateDirectory=outputDir;
designTemplate=[templateDirectory '/template.fsf'];

mkdir(outputDir)

%% Things that need to be looked for in the text file


evTemp=[oldSubjectText '_EV_run*__reg_*']  % format of EV file names
evFormat=[oldSubjectText '_EV_run[1 2 3 4 5 6]__reg_[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30]']

%  provide labels for EVs
%  EVLabels={'outcomeTime', 't=-5'};

% set analysis parameters
brainThresh=.5;

%% Subject loop

cd(rawDatDir)
someFiles=dir('*r01.nii.gz');

%%

if reallyRun==1
    
    % get subject names from the filename of the first run in the data
    % directory
    for i =  1:length(someFiles)
        fn=someFiles(i).name
        suff=strfind(fn, '_r01.nii.gz');
        subName{i}=fn(1:suff-1); % get subject name
        
        IDL=length(subName);
        
        % extract brain
        %anatomicalName= [subName{i}  '_t1' '.nii.gz'];
        exBrainName=[subName{i} '_t1' '_brain' '.nii.gz'];
        %   message=sprintf('bet %s %s -f %f', fullfile(anatDir,anatomicalName), ...
        %       fullfile(anatDir, exBrainName), brainThresh);
        %   [s,w] = unix(message);
        
        
        % find functional data files
        DIR= dir([rawDatDir '/' subName{i} '_r*.nii.gz' ]);
        runFiles={DIR.name};
        
        
        % pause=(5);
        % loop through each data file and do low level data processing and
        % statistics
        for k = 1:length(runFiles)
            
            runFile=fullfile(rawDatDir, runFiles{k});
            
            % parse run filename into relevant info
            newSubName=subName{i}  %****  % need to fix these...
            Run= k
            newRunText=['_r0' num2str(Run)];
            
            %% open design template file
            
            cd(templateDirectory)
            FID=fopen([designTemplate]);                                 % open design template
            allText= fread(FID, 'uint8=>char')';                       %' read design file to character array
            
            %    copyText=allText;
            
            %% Find and replace run information in EV files
            
            newEvName=[evTemp(1:13), num2str(Run), evTemp(15:end-1)];   % create new EV file name beginning based on run
            strLength=length(newEvName) ;                               % get length (for deciding how many characters to write)
            analInd=regexp(allText, evFormat);                          % find start indices of EV names in design template
            
            % replace each instance with EV name that includes correct run
            % information (but old subject name...)
            for kk = 1:length(analInd)
                allText(analInd(kk):analInd(kk)+strLength-1)=newEvName;
            end
            
            
            
            %% Find and replace subject ID (everywhere: EV file names, input filename, output file name)
            
            strLength=length(oldSubjectText);
            newLength=length(newSubName);
            analInd=strfind(allText, oldSubjectText);
            
            % replace with new subject name (stored in newSubName)
            for kk = 1:length(analInd)
                allText=[allText(1:analInd(kk)-1) ...
                    newSubName ...
                    allText(analInd(kk)+strLength:end)];
                analInd=analInd+newLength-strLength;
                
            end
            
            
            %% Find and replace subject number in data file
            strLength=length(oldRunText);
            analInd=strfind(allText, oldRunText);
            
            % replace with new subject name
            for kk = 1:length(analInd)
                allText(analInd(kk):analInd(kk)+strLength-1)=newRunText;
            end
            
            %% Create a new design file
            
            % make unique filename
            
            newDF=[outputDir '/' newSubName num2str(Run) '.fsf'];
            % create file
            newFID=fopen(newDF, 'w+');
            % write the design that
            COUNT=fwrite(newFID, allText);
            
        end
    end
    
end



