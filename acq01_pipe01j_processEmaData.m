% ########## Process EMA Data ##########
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
if exist('GD1','var')
    disp('GD1 already in workspace. Processing subject EMA...')
elseif exist([folderPathExperiment '\' 'groupData.mat'],'file')
    load([folderPathExperiment '\' 'groupData.mat'],'GD1');
end

% Get subject's exported .mat files
cd(folderPathSubjectData)
subjectDataFileList = dir('*.mat');
nFiles = size(subjectDataFileList,1);

% Get subject COR files
cd(folderPathSubjectResults)
corHipFileList = dir('*CorHipRight_pos.sto');
corKneeFileList = dir('*CorKneeRight_pos.sto');
corAnkleFileList = dir('*CorAnkleRight_pos.sto');

% Load subject data table to get peak power
dataTable = readtable([folderPathExperiment '\dataTable.xlsx']);
iSubject = find(contains(dataTable.subject,subjectName));
subjectPeakPower = dataTable.peakPower(iSubject);
subjectMass = dataTable.mass(iSubject);

jointsList = {'hip','knee','ankle'};
nJoints = numel(jointsList);

%% ========== Part 2: Loop through each subject file ==========

% Display progress of script
h = waitbar(0,...
    ['Processing ' subjectName ' EMA Data.']);

% ----- iFiles -----
for iFiles = 1:nFiles
    
    % Update progress
    progress = iFiles/nFiles;
    j = num2str(iFiles);
    n = num2str(nFiles);

    % Compute percentage of completion
    waitbar(progress,h,...
        ['Processing ' subjectName ' EMA Data. File ' j ' of ' n])

    % Set trial name
    trialName = strrep(subjectDataFileList(iFiles).name,'.mat','');

    % ========== Part 2a: Use CoR and MA to calculate EMA ==========
    
    % Get trial struct from COR file
    CorHip = importdata(corHipFileList(iFiles).name,'\t');
    CorKnee = importdata(corKneeFileList(iFiles).name,'\t');
    CorAnkle = importdata(corAnkleFileList(iFiles).name,'\t');

    % Get XY coordinates from COR data
    dataCorHip = [CorHip.data(:,2) CorHip.data(:,3)];
    dataCorKnee = [CorKnee.data(:,2) CorKnee.data(:,3)];
    dataCorAnkle = [CorAnkle.data(:,2) CorAnkle.data(:,3)];

    % Get resultant force variables
    cd(folderPathSubjectResults)
    load([trialName 'workspaceExternalLoads.mat'],'angleClockwiseRadians',...
        'rTan','rRad','lTan','lRad','nFrames','frameRate','reactionForceXyzGlobalRight',...
        'reactionForceXyzGlobalLeft','resultantReactionForceXyzGlobalRight',...
        'resultantReactionForceXyzGlobalLeft','pointXyzGlobalRight','pointXyzGlobalLeft');
    
    % Get resultant force vector in global reference frame
    resultantForceOrigin = pointXyzGlobalRight(1:2,:)';
    resultantForceMagnitude = pointXyzGlobalRight(1:2,:)' + ...
        reactionForceXyzGlobalRight(1:2,:)';

    % Get muscle moment arm lengths
    cd(folderPathSubjectResults)
    hipMomentArms = importdata([trialName ...
        '_MuscleAnalysis_MomentArm_hip_flexion_r.sto'],'\t');
    kneeMomentArms = importdata([trialName ...
        '_MuscleAnalysis_MomentArm_knee_angle_r.sto'],'\t');
    ankleMomentArms = importdata([trialName ...
        '_MuscleAnalysis_MomentArm_ankle_angle_r.sto'],'\t');
    
    % Set muscle moment arms
    rGM = mean(abs(hipMomentArms.data(:,16:18)),2);
    rVL = abs(kneeMomentArms.data(:,40));
    rSOL = abs(ankleMomentArms.data(:,35));
    
    % Calculate perpendicular distance from joint centre to resultant force
    hipR = distancePointLine(dataCorHip,...
        [resultantForceOrigin(1:size(dataCorHip,1),:)...
        resultantForceMagnitude(1:size(dataCorHip,1),:)]);
    
    kneeR = distancePointLine(dataCorKnee,...
        [resultantForceOrigin(1:size(dataCorKnee,1),:)...
        resultantForceMagnitude(1:size(dataCorKnee,1),:)]);
    
    ankleR = distancePointLine(dataCorAnkle,...
        [resultantForceOrigin(1:size(dataCorAnkle,1),:)...
        resultantForceMagnitude(1:size(dataCorAnkle,1),:)]);

    % Calculate EMA
    hipEMA = rGM ./ hipR;
    kneeEMA = rVL ./ kneeR;
    k = kneeEMA > max(rVL / .001);
    kneeEMA(k) = max(rVL / .001);
    ankleEMA = rSOL ./ ankleR;
    
    % ========== Part 2b: Cut data into cycles ==========

    % Find peak locations in angle data
    if contains(trialName,'70')
        targetCadence = 70;
    else
        targetCadence = 120;
    end
    
    % Give 15% buffer for initial peak identification
    buffer = 1.15;
    minPeakDistance = frameRate * (60 / (targetCadence * buffer));
    [pks,locs] = ...
        findpeaks(angleClockwiseRadians,'minpeakdistance',minPeakDistance);
    
    % Create variable to index the number of crank cycles
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
    if contains(trialName,'100')
        targetPower = subjectPeakPower * 0.5;
    else
        targetPower = subjectPeakPower * 0.25;
    end
    
    % Set condition based on trial name
    condition = strrep(trialName(10:end),'_','');
    
    % ----- iJoints -----
    for iJoints = 1:nJoints
        j = 0;
        
        % ----- iLocs -----
        for iLocs = 1:nCycles
            x = 0:1 / (locs(iLocs + 1) - locs(iLocs)):1;
            if strcmp(jointsList(iJoints),'hip')
                v = hipEMA(locs(iLocs):locs(iLocs + 1));
            elseif strcmp(jointsList(iJoints),'knee')
                v = kneeEMA(locs(iLocs):locs(iLocs + 1));
            else
                v = ankleEMA(locs(iLocs):locs(iLocs + 1));
            end
            
            xq = 0:1/100:1;
            vq1 = interp1(x,v,xq,'spline');
            span = 10;            
            
            % Set cut offs for valid cadence and power (\pm 5%)
            lowCut = 0.95;
            highCut = 1.05;
            
            if meanPower(iLocs) > targetPower * lowCut...
                    && meanPower(iLocs) < targetPower * highCut...
                    && meanCadence(iLocs) > targetCadence * lowCut...
                    && meanCadence(iLocs) < targetCadence * highCut
                
                j = j+1;
                EMAcycle = smooth(vq1,span,'rlowess');
                
                if strcmp(jointsList(iJoints),'hip')
                    GD1.(condition).(subjectName).EMA.hip(j,:) = EMAcycle;
                elseif strcmp(jointsList(iJoints),'knee')
                    GD1.(condition).(subjectName).EMA.knee(j,:) = EMAcycle;
                else
                    GD1.(condition).(subjectName).EMA.ankle(j,:) = EMAcycle;
                end
            else
            end
        end
        % ----- iLocs -----
    end
    % ----- iJoints -----
    
    %% ========== Part 3: Save workspace ==========

    % Save moment arm data and EMA results
    cd(folderPathSubjectResults)
    fileName = [trialName 'workspaceEmaAnalysis.mat'];
    save(fileName,'rGM','rVL','rSOL','hipR','kneeR',...
        'ankleR','hipEMA','kneeEMA','ankleEMA')
end
% ----- iFiles -----
    
% Update message
waitbar(1,h,'Saving data structure')

%% ========== Part 4: Save structure ==========
cd(folderPathExperiment)
save('groupData','GD1')

% Update message
waitbar(1,h,[subjectName ' EMA analysis complete.'])
close(h)
%% End