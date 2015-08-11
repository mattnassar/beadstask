%#ok<*SAGROW>
%#ok<*AGROW>
clear classes
clc
currentdir = pwd;
cd '/Users/Carmenere/Documents/ANALYSIS/Notes, Data and Scripts/DATA/BeadsTask/beadstask5 data'
everything = getFileNames('/Users/Carmenere/Documents/ANALYSIS/Notes, Data and Scripts/DATA/BeadsTask/beadstask5 data');

for z = 1:length(everything)
    if isempty(strfind(everything{z},'.'))
        cd (everything{z})
        allfiles = getFileNames(pwd,'.mat');
        for y = 1:length(allfiles)
            data = [];
            if ~isempty(strfind(allfiles{y},'realScannerGame'))
                data = load(allfiles{y});
                data = data.statusData;
            end
            if ~isempty(data)
                setupViewingPhase = [];
                preInfoJitterPhase = [];
                tokenViewingPhase = [];
                postInfoJitterPhase = [];
                freeDrawingPhase = [];
                betNowPhase = [];
                decisionPhase = [];
                tooSlowPhase =[];
                ITIPhase = [];
                leftButtonPressTiming =[];
                rightButtonPressTiming =[];
                for x = 1:length(data)
                    if isfinite(data(x).startTime)
                        keyboard
                        % 3 column regressor: begin time, duration, value
                        startTime = data(x).startTime;
                        setupViewingPhase(end+1,:) = [data(x).preChoiceOn-startTime data(x).preChoiceOff-data(x).preChoiceOn 1];
                        if isfinite(data(x).realInfoOn) % forced info
                            preInfoJitterPhase(end+1,:) = [data(x).realInfoOn-startTime data(x).curr_PostChoiceIs 1];
                            tokenViewingPhase(end+1,:) = [data(x).realInfoOn+data(x).curr_PostChoiceIs-startTime 1 1];
                            postInfoJitterPhase(end+1,:) = [data(x).infoOn-startTime data(x).infoOff-data(x).infoOn 1];
                            if isfinite(data(x).forceChoiceTime)
                                freeDrawingPhase(end+1,:) = [data(x).freeDrawOn-startTime data(x).forceChoiceTime-data(x).freeDrawOn 1];
                                betNowPhase(end+1,:) = [data(x).forceChoiceTime-startTime 0.5 1];
                            elseif isfinite(data(x).tooSlowOnTime) && ~isfinite(data(x).betOn)
                                freeDrawingPhase(end+1,:) = [data(x).freeDrawOn-startTime data(x).tooSlowOnTime-data(x).freeDrawOn 1];
                            else
                                freeDrawingPhase(end+1,:) = [data(x).freeDrawOn-startTime data(x).enterDecTime-data(x).freeDrawOn 1];
                            end
                        end
                        if isfinite(data(x).betOn)
                            if isfinite(data(x).tooSlowOnTime)
                                decisionPhase(end+1,:) = [data(x).betOn-startTime data(x).tooSlowOnTime-data(x).betOn 1];
                            else
                                decisionPhase(end+1,:) = [data(x).betOn-startTime data(x).betchoicetime-data(x).betOn 1];
                            end
                        end
                        if isfinite(data(x).tooSlowOnTime)
                            tooSlowPhase(end+1,:) =[data(x).tooSlowOnTime-startTime data(x).tooSlowOffTime-data(x).tooSlowOnTime 1];
                            ITIPhase(end+1,:) = [data(x).tooSlowOffTime-startTime data(x).curr_InterTrialIs 1];
                        else
                            ITIPhase(end+1,:) = [data(x).betchoicetime-startTime data(x).curr_InterTrialIs 1];
                        end
                        if ~isempty(data(x).extraDrawTime) && sum(isfinite(data(x).extraDrawTime))>0
                            if data(x).infoButtonSide == 0 %info button was on the left
                                leftButtonPressTiming =[leftButtonPressTiming data(x).extraDrawTime-startTime];
                                if ~isfinite(data(x).forceChoiceTime)
                                    rightButtonPressTiming = [rightButtonPressTiming data(x).enterDecTime-startTime];
                                end
                            else
                                rightButtonPressTiming =[rightButtonPressTiming data(x).extraDrawTime-startTime];
                                if ~isfinite(data(x).forceChoiceTime)
                                    leftButtonPressTiming = [leftButtonPressTiming data(x).enterDecTime-startTime];
                                end
                            end
                        end
                        if isfinite(data(x).betchoicetime)
                            if data(x).curr_choice == 1 %subject chose left option
                                leftButtonPressTiming = [leftButtonPressTiming data(x).betchoicetime-startTime];
                            else
                                rightButtonPressTiming = [rightButtonPressTiming data(x).betchoicetime-startTime];
                            end
                        end
                    end
                end
                buttonPressTiming = [leftButtonPressTiming' ; rightButtonPressTiming'];
                buttonPressTiming = sort(buttonPressTiming);
                buttonPressTiming = [buttonPressTiming, ones(size(buttonPressTiming)).*.1 ones(size(buttonPressTiming))];
                leftButtonPressTiming = [leftButtonPressTiming' ones(size(leftButtonPressTiming')).*.1 ones(size(leftButtonPressTiming')).*(-1)];
                rightButtonPressTiming = [rightButtonPressTiming' ones(size(rightButtonPressTiming')).*.1 ones(size(rightButtonPressTiming'))];
                pressRightButton = sortrows([leftButtonPressTiming; rightButtonPressTiming]);
                pressRightButton(:,3) = zscore(pressRightButton(:,3));
                try
                    cd 'EVs'
                catch
                    mkdir 'EVs'
                    cd 'EVs'
                end
                dlmwrite([allfiles{y}(1:27) 'setupViewingPhase'],setupViewingPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'preInfoJitterPhase'],preInfoJitterPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'tokenViewingPhase'],tokenViewingPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'postInfoJitterPhase'],postInfoJitterPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'freeDrawingPhase'],freeDrawingPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'betNowPhase'],betNowPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'decisionPhase'],decisionPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'tooSlowPhase'],tooSlowPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'ITIPhase'],ITIPhase,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'buttonPressTiming'],buttonPressTiming,'delimiter','\t')
                dlmwrite([allfiles{y}(1:27) 'pressRightButton'],pressRightButton,'delimiter','\t')
                cd ..
            end
        end
        cd ..
    end
end

cd (currentdir)