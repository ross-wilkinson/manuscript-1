% ########## Process Kinetic Data ##########
% 
% Description
% ---------------
% This script processes the inverse dynamics output from OpenSim 3.3. Data
% is first cut into crank cycles, interpolated to 101 points and then valid
% data i.e. within target cadence and power thresholds is placed into a
% group data structure (GD1).
%
% Pre-requisites
% --------------- 
% Nil - besides having ID .sto file outputs from OpenSim 3.3
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
    disp('GD1 already in workspace. Processing subject kinetics...')
elseif exist([folderPathExperiment '\' 'groupData.mat'],'file')
    load([folderPathExperiment '\' 'groupData.mat'],'GD1');
end
    
% Find subject's exported .mat files
cd(folderPathSubjectData)
subjectDataFileList = dir('*.mat');
nFiles = size(subjectDataFileList,1);

% Find subject id files
cd(folderPathSubjectResults)
subjectIdFileList = dir('*inverseDynamics.sto');

% Load subject data table to find peak power
dataTable = readtable([folderPathExperiment '\dataTable.xlsx']);
iSubject = find(contains(dataTable.subject,subjectName));
subjectPeakPower = dataTable.peakPower(iSubject);
subjectMass = dataTable.mass(iSubject);

%% ========== Part 2: Loop through each subject file ==========

% Display progress of script
h = waitbar(0,...
    ['Processing ' subjectName ' Kinetic Data.']);

% ----- iFiles -----
for iFiles = 1:nFiles
    
    % Update progress
    progress = iFiles/nFiles;
    j = num2str(iFiles);
    n = num2str(nFiles);

    % Compute percentage of completion
    waitbar(progress,h,...
        ['Processing ' subjectName ' Kinetic Data. File ' j ' of ' n])
    
    % Load data file
    cd(folderPathSubjectData)
    data = load(subjectDataFileList(iFiles).name);
    trialName = fieldnames(data);
    
    % ========== Part 2a: Cut data into cycles ==========
    
    % Load id results file and some of the force data workspace variables
    cd(folderPathSubjectResults)
    dataInverseDynamics = importdata(subjectIdFileList(iFiles).name,'\t');
    load([trialName{1} 'workspaceExternalLoads.mat'],'angleClockwiseRadians',...
        'rTan','rRad','lTan','lRad','nFrames','frameRate',...
        'reactionForceXyzGlobalRight','reactionForceXyzGlobalLeft',...
        'resultantReactionForceXyzGlobalRight',...
        'resultantReactionForceXyzGlobalLeft','pointXyzGlobalRight',...
        'pointXyzGlobalLeft');
    
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
    
    % Create joint moment data and label variables
    jointMomentData = dataInverseDynamics.data;
    jointMomentLabels = categorical(dataInverseDynamics.colheaders);
    
    % Set condition based on file name.
    condition = strrep(trialName{1}(10:end),'_','');

    % Gather force data for loop within loop
    forceData = ...
        vertcat(rTan',rRad',lTan',lRad',resultantReactionForceXyzGlobalRight,...
    resultantReactionForceXyzGlobalLeft,reactionForceXyzGlobalRight,...
    reactionForceXyzGlobalLeft,pointXyzGlobalRight,pointXyzGlobalLeft);
    
    forces = {'rTan','rRad','lTan','lRad',...
        'resultantReactionForceXyzGlobalRight',...
        'resultantReactionForceXyzGlobalLeft','reactionForceGlobalRightX',...
        'reactionForceGLobalRightY','reactionForceGLobalRightZ',...
        'reactionForceGlobalLeftX','reactionForceGLobalLeftY',...
        'reactionForceGLobalLeftZ','pointRightX','pointRightY','pointRightZ',...
        'pointLeftX','pointLeftY','pointLeftZ'};
    
    % ========== Part 2b: Interpolate each cycle ==========
    
    % ----- iCols -----
    for iCols = 2:size(jointMomentData,2)
        
        % Set data
        colData = jointMomentData(:,iCols);
        
        % Valid data counter
        j = 0;
        
        % ----- iLocs -----
        for iLocs = 1:nCycles
            x = 0:1 / (locs(iLocs + 1) - locs(iLocs)):1;
            v = colData(locs(iLocs):locs(iLocs + 1))';
            xq = 0:1/100:1;
            vq1 = interp1(x,v,xq,'spline');
            if cadence == 120
                span = 19;
            else
                span = 11;
            end
            jointMomentInterpolated = smooth(vq1,span,'sgolay');
            
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
                GD1.(condition).(subjectName).jointMoment. ...
                    (char(jointMomentLabels(iCols)))(j,:) = ...
                    jointMomentInterpolated;
                ExtNmKg = [char(jointMomentLabels(iCols)) 'PerKg'];
                GD1.(condition).(subjectName).jointMoment.(ExtNmKg)(j,:) = ...
                    jointMomentInterpolated / subjectMass;
                
                % Skip pelvis force variables
                switch iCols
                    case {5,6,7} 
                    otherwise
                        jointVelocityLabel = ...
                            strrep(string(jointMomentLabels(iCols)),...
                            '_moment','');
                        jointVelocity = ...
                            GD1.(condition).(subjectName).jointVelocity. ...
                            (jointVelocityLabel)(j,:) / (180 / pi);
                        GD1.(condition).(subjectName).jointPower. ...
                            (jointVelocityLabel)(j,:) = ...
                            jointMomentInterpolated(:)' .* jointVelocity;
                        jointPower = ...
                            GD1.(condition).(subjectName).jointPower. ...
                            (jointVelocityLabel)(j,:);
                        ExtWKg = [char(jointVelocityLabel) 'PerKg'];
                        GD1.(condition).(subjectName).jointPower. ...
                            (ExtWKg)(j,:) = jointPower / subjectMass;
                        
                        spacing = 60 / 100 / cadence;
                        % Calculate -ve and +ve joint power
                        pwr = GD1.(condition).(subjectName).jointPower. ...
                            (ExtWKg)(j,:);
                        % +ve velocity = flexion
                        vel = jointVelocity;
                        % +ve moment = flexion
                        mom = GD1.(condition).(subjectName).jointMoment. ...
                            (ExtNmKg)(j,:);
                        negativeExtensionLabel = ...
                            [char(jointVelocityLabel) 'NegExt'];
                        negativeFlexionLabel = ...
                            [char(jointVelocityLabel) 'NegFlex'];
                        positiveExtensionLabel = ...
                            [char(jointVelocityLabel) 'PosExt'];
                        positiveFlexionLabel = ...
                            [char(jointVelocityLabel) 'PosFlex'];
                        % pos Flexor power
                        GD1.(condition).(subjectName).jointPower. ...
                            (positiveFlexionLabel)(j) = ...
                            trapz(pwr(vel > 0 & mom > 0)) * spacing; 
                        % pos Extensor power
                        GD1.(condition).(subjectName).jointPower. ...
                            (positiveExtensionLabel)(j) = ...
                            trapz(pwr(vel < 0 & mom < 0)) * spacing;
                        % neg Flex power (Flexor moment + Ext velocity)
                        GD1.(condition).(subjectName).jointPower. ...
                            (negativeFlexionLabel)(j) = ...
                            trapz(pwr(vel < 0 & mom > 0)) * spacing; 
                        % neg Ext power (Extensor moment + Flex velocity) 
                        GD1.(condition).(subjectName).jointPower. ...
                            (negativeExtensionLabel)(j) = ...
                            trapz(pwr(vel > 0 & mom < 0)) * spacing;       
                         
                        % Normalise joint work to body mass
                        jointWork = trapz(jointPower) * spacing;
                        GD1.(condition).(subjectName).jointWork. ...
                            (jointVelocityLabel)(j) = jointWork;
                        ExtWorkKg = [char(jointVelocityLabel) 'PerKg'];
                        GD1.(condition).(subjectName).jointWork. ...
                            (ExtWorkKg)(j) = jointWork / subjectMass;
                        jointWorkKg = GD1.(condition).(subjectName). ...
                            jointWork.(ExtWorkKg)(j);
                        % Cumulative work in a minute
                        ExtWorkKgMin = [char(jointVelocityLabel) 'PerKgMin'];
                        GD1.(condition).(subjectName).jointWork. ...
                            (ExtWorkKgMin)(j) = jointWorkKg * 60;   
                end
            else
            end
            
            % Cut force data into crank cycles
            % ----- iForces -----
            for iForces = 1:size(forceData,1)
                f = forceData(iForces,locs(iLocs):locs(iLocs + 1));
                fq1 = interp1(x,f,xq,'spline');
                
                if meanPower(iLocs) > targetPower * lowCut...
                        && meanPower(iLocs) < targetPower * highCut...
                        && meanCadence(iLocs) > cadence * lowCut...
                        && meanCadence(iLocs) < cadence * highCut
                    GD1.(condition).(subjectName).crankForce. ...
                        (forces{iForces})(j,:) = fq1;
                else
                end
                
            end
            % ----- iForces -----
        end
        % ----- iLocs -----
    end
    % ----- iCols -----
end
% ----- iFiles -----

%% ========== Part 3: Save structure & plots ==========

% Update message
waitbar(1,h,'Saving data structure')

% Save group data structure
save([folderPathExperiment '\' 'groupData'],'GD1')

% Update message
waitbar(1,h,[subjectName ' Kinetics data processing complete.'])
close(h)

%% End