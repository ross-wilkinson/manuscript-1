% ########## Process Kinematic Data ##########
% 
% Description
% ---------------
% This script processes the inverse kinematic output from OpenSim 3.3. Data
% is first cut into crank cycles, interpolated to 101 points and then valid
% data i.e. within target cadence and power thresholds is placed into a
% group data structure (GD1).
%
% Pre-requisites
% --------------- 
% Nil - besides having IK .mot file outputs from OpenSim 3.3
% 
% Acknowledgement
% ----------------
% Written by Ross Wilkinson, The University of Queensland,
% (ross.wilkinson_at_uqconnect.edu.au)

%% Initialise
clearvars -except GD1;clc;close all

%% ========== Part 1: Get and set subject data and folders ==========

% Input subject number
subjectNo = input('Subject Number? ','s');

% Create subject name
subjectName = ['subject' subjectNo];

% Set filepath for subject setup and results files
folderPathExperiment = 'D:\exp01';
folderPathSubjectSetup = [folderPathExperiment '\' subjectName '\setup'];
folderPathSubjectResults = [folderPathExperiment '\' subjectName '\results'];
folderPathSubjectData = [folderPathExperiment '\' subjectName '\data'];

% Load group data structure if it exists
if exist('GD1','var') == 1
    disp('GD1 already in workspace. Processing subject kinematics...')
elseif exist([folderPathExperiment '\' 'groupData.mat'],'file')
    load([folderPathExperiment '\' 'groupData.mat'],'GD1');
end

% Find subject's exported .mat files
cd(folderPathSubjectData)
subjectDataFileList = dir('*.mat');
nFiles = size(subjectDataFileList,1);

% Find subject ik files
cd(folderPathSubjectResults)
subjectIkFileList = dir('*inverseKinematics.mot');

% Load subject data table to find peak power
dataTable = readtable([folderPathExperiment '\dataTable.xlsx']);
iSubject = find(contains(dataTable.subject,subjectName));
subjectPeakPower = dataTable.peakPower(iSubject);

% Open figure to plot power comparison
figure('units','normalized','position',[0.1 0.1 .8 .8]);
subplot(2,3,nFiles)

%% ========== Part 2: Loop through each subject file ==========

% Display progress of script
h = waitbar(0,...
    ['Processing ' subjectName ' Kinematic Data.']);

% ----- iFiles -----
for iFiles = 1:nFiles
    
    % Update progress
    progress = iFiles/nFiles;
    j = num2str(iFiles);
    n = num2str(nFiles);

    % Compute percentage of completion
    waitbar(progress,h,...
        ['Processing ' subjectName ' Kinematic Data. File ' j ' of ' n])

    % Load data file
    cd(folderPathSubjectData)
    data = load(subjectDataFileList(iFiles).name);
    trialName = fieldnames(data);

    % ========== Part 2a: Cut data into cycles ==========

    % Load ik results file and some of the force data workspace variables
    cd(folderPathSubjectResults)
    dataInverseKinematics = importdata(subjectIkFileList(iFiles).name,'\t');
    load([trialName{1} 'workspaceExternalLoads.mat'],...
        'angleClockwiseRadians','rTan','nFrames','frameRate');

    % Find peaks in crank angle signal based off cadence and frame rate
    if contains(trialName{1},'70')
        cadence = 70;
    else
        cadence = 120;
    end

    % Give 15% buffer for initial peak identification
    buffer = 1.15;
    minPeakDistance = frameRate * (60 / (cadence * buffer));
    [pks,locs] = ...
        findpeaks(angleClockwiseRadians,'minpeakdistance',minPeakDistance);

    % Calculate average power for each pedal cycle on right crank.
    nCycles = length(locs) - 1;
    [meanCadence,meanVelocity,meanTorque,meanPower] = deal(zeros(1,nCycles));
    crankLength = 0.175;
    
    % ----- iLocs -----
    for iLocs = 1:nCycles
        meanCadence(iLocs) = 60 * frameRate / (locs(iLocs + 1) - locs(iLocs));
        meanVelocity(iLocs) = meanCadence(iLocs) / 60 * 2 * pi;
        meanTorque(iLocs) = ...
            mean(rTan(locs(iLocs):locs(iLocs + 1))) * crankLength;
        meanPower(iLocs) = meanTorque(iLocs) * meanVelocity(iLocs);
    end
    % ----- iLocs -----

    % Set target power based on file name. Halve as only for right side
    if contains(trialName{1},'100')
        targetPower = subjectPeakPower * 0.5;
    else
        targetPower = subjectPeakPower * 0.25;
    end

    % Create joint angle data and label variables
    jointAngleData = dataInverseKinematics.data;
    jointAngleLabels = categorical(dataInverseKinematics.colheaders);

    % Set condition based on file name.
    condition = strrep(trialName{1}(10:end),'_','');

    % ========== Part 2b: Interpolate each cycle ==========

    % Set step size for finding joint velocity using 'diff' function
    stepTime = 60 / cadence / 100; % step size

    % ----- iCols -----
    for iCols = 2:size(jointAngleData,2)

        % Set data
        colData = jointAngleData(:,iCols);

        % Valid data counter
        j = 0;

        % ----- iLocs -----
        for iLocs = 1:nCycles
            x = 0:1 / (locs(iLocs + 1) - locs(iLocs)):1;
            v = colData(locs(iLocs):locs(iLocs + 1))';
            xq = 0:1 / 100:1;
            vq1 = interp1(x,v,xq,'spline');
            if cadence == 120
                span = 19;
            else
                span = 11;
            end
            jointAngleInterpolated = smooth(vq1,span,'sgolay');

            % Set cut offs for valid cadence and power (\pm 5%)
            lowCut = 0.95;
            highCut = 1.05;

            % ========== Part 2c: Add valid data to structure ==========

            % Use cut offs to assign data to structure
            if meanPower(iLocs) > targetPower * lowCut...                    
                    && meanPower(iLocs) < targetPower * highCut...
                    && meanCadence(iLocs) > cadence * lowCut...
                    && meanCadence(iLocs) < cadence * highCut
                
                j = j+1;
                GD1.(condition).(subjectName).jointAngle. ...
                    (string(jointAngleLabels(iCols)))(j,:) = ...
                    jointAngleInterpolated;
                jointVelRaw = diff(jointAngleInterpolated)/stepTime;
                jointVelSmooth = smooth(jointVelRaw,span,'sgolay')';
                nInt = length(jointVelSmooth);
                jointVelInterpolated = ...
                    interp1(1:nInt,jointVelSmooth,1:(nInt-1)/nInt:nInt);
                GD1.(condition).(subjectName).jointVelocity. ...
                    (string(jointAngleLabels(iCols)))(j,:) = ...
                    jointVelInterpolated;
            else
            end
        end
        % ----- iLocs -----
    end
    % ----- iCols -----

    j = 0;
    % ----- iLocs -----
    for iLocs = 1:nCycles
        if meanPower(iLocs) > targetPower * lowCut...                    
                    && meanPower(iLocs) < targetPower * highCut...
                    && meanCadence(iLocs) > cadence * lowCut...
                    && meanCadence(iLocs) < cadence * highCut
            j = j+1;
            validMeanPowerData(j) = meanPower(iLocs);
            validMeanCadenceData(j) = meanCadence(iLocs);
        else
        end
    end
    % ----- iLocs -----
    
    % Save power output and cadence values for each valid cycle
    GD1.(condition).(subjectName).crankPower.powerTarget = targetPower;

    GD1.(condition).(subjectName).crankPower.powerCycle = validMeanPowerData;

    GD1.(condition).(subjectName).crankPower.powerMean = ...
        mean(validMeanPowerData);

    GD1.(condition).(subjectName).crankPower.powerStd = ...
        std(validMeanPowerData);

    GD1.(condition).(subjectName).crankPower.powerPercDiff = ...
       (GD1.(condition).(subjectName).crankPower.powerMean - ...
       GD1.(condition).(subjectName).crankPower.powerTarget) / ...
       GD1.(condition).(subjectName).crankPower.powerTarget * 100;

    GD1.(condition).(subjectName).cadence.cadenceCycle = validMeanCadenceData;

    GD1.(condition).(subjectName).cadence.cadenceMean = ...
        mean(validMeanCadenceData);

    GD1.(condition).(subjectName).cadence.cadenceStd = ...
        std(validMeanCadenceData);

    GD1.(condition).(subjectName).cadence.cadencePercDiff = ...
       (GD1.(condition).(subjectName).cadence.cadenceMean-cadence)...
       / cadence * 100;

    % Create sub plots of average power against target power threshold for
    % each trial
    ax = subplot(2,3,iFiles);
    hold on
    title(ax,trialName{1},'Interpreter','none');
    ylabel('power (W)')
    ylim([0 subjectPeakPower * 0.6])
    yyaxis right
    hold on
    ylabel('cadence (rpm)')
    ylim([50 150])
    
    if exist('validMeanPowerData','var')
        yyaxis left
        bar(ax,validMeanPowerData);
        line([1 length(validMeanPowerData)],...
            [targetPower targetPower],'Color','blue','Linestyle','--')
        yyaxis right
        plot(ax,validMeanCadenceData,'r*');
        line([1 length(validMeanPowerData)],...
            [cadence cadence],'Color','red','Linestyle','--')
        clear validMeanPowerData validMeanCadenceData
    else
    end
end  
% ----- iFiles -----

%% ========== Part 3: Save structure & plots ==========

% Update message
waitbar(1,h,'Saving data structure')

% Save group data structure
save([folderPathExperiment '\' 'groupData'],'GD1')

% Update message
waitbar(1,h,[subjectName ' Kinematics data processing complete.'])

% Save subject power plots
savefig([folderPathSubjectResults '\' subjectName 'plotCheckPower']);
close all

%% ========== Part 4: Visualise valid data ==========

% Update message
waitbar(1,h,[subjectName ' Plotting kinematic data.'])

% Plot power and cadence comparison for subject
figure
x = 1:4;
y = [GD1.seated5070.(subjectName).crankPower.powerPercDiff ...
    GD1.standing5070.(subjectName).crankPower.powerPercDiff...
    GD1.seated50120.(subjectName).crankPower.powerPercDiff...
    GD1.standing50120.(subjectName).crankPower.powerPercDiff];
bar(x,y);
hold on;
line([0 5],[0 0],'Color','r')
xlabel('Condition')
ylabel('% away from target (watts)')
xticklabels({'S70','NS70','S120','NS120'})
title(['Power comparison ' subjectName ' 50%'])

savefig([folderPathSubjectResults '\' subjectName 'plotComparePowerVsTarget']);
close all

figure
y = [GD1.seated5070.(subjectName).cadence.cadencePercDiff...
    GD1.standing5070.(subjectName).cadence.cadencePercDiff...
    GD1.seated50120.(subjectName).cadence.cadencePercDiff...
    GD1.standing50120.(subjectName).cadence.cadencePercDiff];
bar(x,y);
hold on;
%line([0 5],[0 0],'Color','red')
xlabel('Condition')
ylabel('% away from target (rpm)')
xticklabels({'S70','NS70','S120','NS120'})
title(['Cadence comparison ' subjectName ' 50%'])

savefig([folderPathSubjectResults '\' subjectName ...
    'plotCompareCadenceVsTarget']);

close(h)
close all
%% End