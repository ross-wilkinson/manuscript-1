% ########## Write External Loads File ##########
% 
% Description
% ----------------
% This script writes a subject specific external loads file to be used
% within the Inverse Dynamics setup file needed to run the Inverse Dynamics
% Tool in OpenSim 3.3.
%
% Pre-requisites
% ---------------- 
% To execute this script you will need to:
%   1. Download the MATLAB-OPENSIM INTERFACE TOOLBOX which can be accessed
%   here: https://simtk.org/projects/matlab_tools
%   2. Create a generic structure file which has many of the structure
%   fields pre-filled with default settings. This can be done by running
%   the 'writeSetupInverseKinematicsXml.m' script, which is available within
%   my other GitHub repo: https://github.com/uqrwilk/opensimModelling
%
% Acknowledgement
% ----------------
% This script is merely an adaptation of Glen Lichtwark's (The University
% of Queensland) MatlabOpensimPipelineTools 'grf2xml.m' function. The
% functions available within the above mentioned toolbox are far more
% efficient and robust than my scripts, however I made these scripts as a
% learning process for myself.
% 
% Written by Ross Wilkinson, The University of Queensland,
% (ross.wilkinson_at_uqconnect.edu.au)

%% Initialise
clear;clc;close all
import org.opensim.modeling.*

%% ========== Part 1: Get and set subject data and folders ==========

% Load structure of generic External Loads file
load('D:\opensim\structureExternalLoads.mat')

% Input subject number
subjectNo = input('Subject Number? ','s');

% Create subject name
subjectName = ['subject' subjectNo];

% Set filepath for subject setup and results files
folderPathExperiment = 'D:\exp01';
folderPathSubjectSetup = [folderPathExperiment '\' subjectName '\setup'];
folderPathSubjectResults = [folderPathExperiment '\' subjectName '\results'];
folderPathSubjectData = [folderPathExperiment '\' subjectName '\data'];

% Find subject's exported .mat files
cd(folderPathSubjectData)
subjectFileList = dir('*.mat');
nFiles = size(subjectFileList,1);

% Find side that crank angle was recorded from in data table
dataTable = readtable([folderPathExperiment '\dataTable.xlsx']);
subjectIndex = find(contains(dataTable.subject,subjectName));
crankSide = dataTable.crankSide{subjectIndex};
samplingFactor = dataTable.samplingFactor{subjectIndex};

% Loop: Process crank force data
disp('Writing External Loads files...')

%% ========== Part 2: Loop through each subject file ==========

% Display progress of script
h = waitbar(0,...
    ['Processing ' subjectName ' external load files.']);

% ----- iFiles -----
for iFiles = 1:nFiles
    
    % Update progress
    progress = iFiles/nFiles;
    j = num2str(iFiles);
    n = num2str(nFiles);
    
    % Fill percentage of waitbar
    waitbar(progress,h,...
        ['Processing ' subjectName ' external load files. File ' j ' of ' n])
    
    % Load file data
    cd(folderPathSubjectData)
    data = load(subjectFileList(iFiles).name);
    trialName = fieldnames(data);
    
    % Set marker data, marker labels and analog labels to vars
    markerData = data.(trialName{1}).Trajectories.Labeled.Data(:,1:3,:);
    markerLabels = categorical(data.(trialName{1}).Trajectories.Labeled.Labels);
    analogLabels = categorical(data.(trialName{1}).Analog.Labels);

    % Re-shape marker_data for writing .trc file. One marker per page.
    markerData = permute(markerData,[2 3 1]);

    % Set specific marker and analog vars
    B2 = markerData(:,:,markerLabels == 'B2');
    B3 = markerData(:,:,markerLabels == 'B3');
    rtoe = markerData(:,:,markerLabels == 'rtoe');
    rmt5 = markerData(:,:,markerLabels == 'rmt5');
    rcal = markerData(:,:,markerLabels == 'rcal');
    ltoe = markerData(:,:,markerLabels == 'ltoe');
    lmt5 = markerData(:,:,markerLabels == 'lmt5');
    lcal = markerData(:,:,markerLabels == 'lcal');
    
    % ========== Part 2a: Process force and angle data ==========
    
    % ----- samplingFactor -----
    switch samplingFactor
        case 'yes' % Take all data
            rTan = data.(trialName{1}).Analog.Data(analogLabels == 'R.Tan',:);
            rRad = data.(trialName{1}).Analog.Data(analogLabels == 'R.Rad',:);
            lTan = data.(trialName{1}).Analog.Data(analogLabels == 'L.Tan',:);
            lRad = data.(trialName{1}).Analog.Data(analogLabels == 'L.Rad',:);
            angle = data.(trialName{1}).Analog.Data(analogLabels == 'Angle',:);          
            % ----- subjectNo -----
            switch subjectNo
                case '08'
                    lTan = data.(trialName{1}).Analog.Data(analogLabels == 'L.Rad',:);
                    lRad = data.(trialName{1}).Analog.Data(analogLabels == 'L.Tan',:);
            end
            % ----- subjectNo -----   
            % Interpolate data to match frame rate
            samplingFactorNew = data.(trialName{1}).Analog.SamplingFactor;
            nSamples = data.(trialName{1}).Analog.NrOfSamples;
            rTan = interp1(rTan,1:samplingFactorNew:nSamples);
            rRad = interp1(rRad,1:samplingFactorNew:nSamples);
            lTan = interp1(lTan,1:samplingFactorNew:nSamples);
            lRad = interp1(lRad,1:samplingFactorNew:nSamples);
            angle = interp1(angle,1:samplingFactorNew:nSamples);
        case 'no' % Only take data from inital number of frames
            rTan = data.(trialName{1}).Analog.Data(analogLabels == ...
                'R.Tan',1:data.(trialName{1}).Frames);
            rRad = data.(trialName{1}).Analog.Data(analogLabels == ...
                'R.Rad',1:data.(trialName{1}).Frames);
            lTan = data.(trialName{1}).Analog.Data(analogLabels == ...
                'L.Tan',1:data.(trialName{1}).Frames);
            lRad = data.(trialName{1}).Analog.Data(analogLabels == ...
                'L.Rad',1:data.(trialName{1}).Frames);
            angle = data.(trialName{1}).Analog.Data(analogLabels == ...
                'Angle',1:data.(trialName{1}).Frames);
            % ----- subjectNo -----
            switch subjectNo
                case '04'
                    angle = data.(trialName{1}).Analog.Data...
                        (37,1:data.(trialName{1}).Frames);
            end
            % ----- subjectNo -----
    end
    % ----- samplingFactor -----

    % Remove offset in force signals
    
    % ----- subjectNo -----
    switch subjectNo
        case {'06','09','20'} 
            offsetTanRight = 2.5;
            offsetRadRight = 2.435;
            offsetTanLeft = 2.5;
            offsetRadLeft = 2.4662;
        case {'04','08','12','23'}
            offsetTanRight = 2.5;
            offsetRadRight = 2.34;
            offsetTanLeft = 2.5;
            offsetRadLeft = 2.435;
        otherwise
            offsetTanRight = 2.5;
            offsetRadRight = 2.5;
            offsetTanLeft = 2.5;
            offsetRadLeft = 2.5;
    end
    % ----- subjectNo -----
    
    rTan = rTan-offsetTanRight;
    rRad = rRad-offsetRadRight;
    lTan = lTan-offsetTanLeft;
    lRad = lRad-offsetRadLeft;

    % Convert signal voltage to SI Units (Newtons).
    crankLength = 0.175; % 175 mm
    mvToForceTan = 200/crankLength; % 5mV/N. 0.1725m crank length
    mvToForceRad = 1000; % 1mV/N
    mvToDeg = 72; % 72 degrees/V
    
    rTan = rTan*mvToForceTan;
    rRad = rRad*mvToForceRad;
    lTan = lTan*mvToForceTan;
    lRad = lRad*mvToForceRad;
    angle = angle*mvToDeg;

    % Check which side crank angle was recorded from. SWITCH if necessary
    % ----- crankSide -----
    switch crankSide
        case 'Left'
            angleOpp = angle - 180;
            angleOpp(angleOpp<0) = angleOpp(angleOpp<0) + 360;
            angle = angleOpp;
        otherwise
    end
    % ----- crankSide -----

    % ========== Part 2b: Visual inspection of crank signals ==========
    
    % Plot signals before filtering and angle calculations in subplots
    
    close
    figure
    ax1 = subplot(3,2,1);
    plot(ax1,lTan);
    title(ax1,'ltan');
    hold on;

    ax2 = subplot(3,2,2);
    plot(ax2,rTan);
    title(ax2,'rtan');
    hold on;

    ax3 = subplot(3,2,3);
    plot(ax3,lRad);
    title(ax3,'lrad');
    hold on;

    ax4 = subplot(3,2,4);
    plot(ax4,rRad);
    title(ax4,'rrad');
    hold on;

    ax5 = subplot(3,2,[5 6]);
    plot(ax5,angle);
    title(ax5,'angle');
    hold on;

    % ========== Part 2c: Coordinate rotation ==========

    % Smooth noise out of data. Set span based on frame rate.
    frameRate = data.(trialName{1}).FrameRate;
    span = frameRate/10;
    filter = 'sgolay';

    rTan = smooth(rTan,span,filter);
    rRad = smooth(rRad,span,filter);
    lTan = smooth(lTan,span,filter);
    lRad = smooth(lRad,span,filter);

    % Calculate z-axis unit vector
    % Pre-allocate point and force arrays
    [pointXyzGlobalRight, pointXyzGlobalLeft] = deal(zeros(3,length(rTan)));
    [forceXyzGlobalRight, forceXyzGlobalLeft] = ...
        deal(zeros(3,length(rTan)));
    [angleClockwiseRadians, angleCounterClockwiseRadians, ...
        angleGlobalRightRadians, angleGlobalLeftRadians] = deal(zeros(1,length(rTan)));
        
    % Displacement from B3 to B2 is new y-axis
    displacement = B2(:,1)' - B3(:,1)';
    % New y-axis in bicycle reference frame. Divide by Euclidean length of the
    %vector.
    yAxisBicycle = displacement / norm(displacement);
    % Use global z-axis to calculate x-axis in bicycle reference frame. 
    zAxisBicycle = [0 0 1];
    % Take cross product of y(bicycle) and z(global) to find new x-axis in
    % bicycle reference frame.
    xAxisBicycle = cross(yAxisBicycle,zAxisBicycle); 
    % Now use x and y-axis in bicycle reference frame to calculate new z-axis.
    zAxisBicycle = cross(xAxisBicycle,yAxisBicycle);
    % Set bicycle unit vectors in rotation matrix
    rotationMatrix3d = [xAxisBicycle;yAxisBicycle;zAxisBicycle];

    % Rotate marker data
    markerDataBicycle = zeros(size(markerData));
    % ----- i -----
    for i = 1:length(markerLabels)
        markerDataBicycle(:,:,i) = rotationMatrix3d * markerData(:,:,i);
    end
    % ----- i -----  
    % Loop through calculations for each data point
    % ----- iData -----
    for iData = 1:length(rTan)
        % Calculate virtual marker positions
        % Measure displacement from calc to toe on each foot
        dispCalcToToeRight = rtoe(:,iData) - rcal(:,iData);
        dispCalcToToeLeft = ltoe(:,iData) - lcal(:,iData);
        % Set virtual marker pos as toe minus 1/3 of the disp. to calc in
        % global ref. frame
        virtualMarkerRight = rtoe(:,iData) - dispCalcToToeRight * (1/3);
        virtualMarkerLeft = ltoe(:,iData) - dispCalcToToeLeft * (1/3);
        % Use mt5 for Z coordinate
        virtualMarkerRight(3) = rmt5(3,iData);
        virtualMarkerLeft(3) = lmt5(3,iData);
        % Calculate crank angle using z-unit & virtual marker
        % Set inputs for angle3Points.
        p1 = [virtualMarkerLeft(1) + zAxisBicycle(1) , virtualMarkerLeft(3) + ...
            zAxisBicycle(3)];
        p2 = [virtualMarkerLeft(1) virtualMarkerLeft(3)];
        p3 = [virtualMarkerRight(1) virtualMarkerRight(3)];
        % Calculate crank angle clockwise from TDC
        angleClockwiseRadians(iData) = 2 * pi - angle3Points(p1,p2,p3);

        % Convert crank angle from clockwise to counter-clockwise
        angleCounterClockwiseRadians(iData) = 2 * pi - angleClockwiseRadians(iData);

        % Convert crank angle to global (x/y) reference
        % frame. Rotated clockwise pi/2 (90 deg) from global.
        angleGlobalRightRadians(iData) = angleCounterClockwiseRadians(iData) + (pi / 2);

        % Use global_angle_r to create left crank angle.
        % Left side of crank always ahead of other by pi (180 deg).
        angleGlobalLeftRadians(iData) = angleGlobalRightRadians(iData) + pi;

        % Calculate global force vector direction using crank angle

        % Create rotation matrix based on crank angle
        thetaRight = angleGlobalRightRadians(iData);
        thetaLeft = angleGlobalLeftRadians(iData);
        rotationMatrix2dRight = ...
            [cos(thetaRight) -sin(thetaRight); sin(thetaRight) cos(thetaRight)];
        rotationMatrix2dLeft = ...
            [cos(thetaLeft) -sin(thetaLeft); sin(thetaLeft) cos(thetaLeft)];
        % Rotate crank forces from crank reference frame to bicycle reference
        % frame
        % Positive x in crank reference frame is equal to -ve radial force
        forceXaxisCrankRight = -rRad(iData); 
        forceXaxisCrankLeft = -lRad(iData);
        % Positive y in crank reference frame is equal to -ve tangential
        % force
        forceYaxisCrankRight = -rTan(iData);
        forceYaxisCrankLeft = -lTan(iData);       
        % Multiply by 2D rotation matrix
        forceXyzBicycleRight = rotationMatrix2dRight * ...
            [forceXaxisCrankRight;forceYaxisCrankRight];
        forceXyzBicycleLeft = rotationMatrix2dLeft * ...
            [forceXaxisCrankLeft;forceYaxisCrankLeft];
        % Set z coordinate at zero. No z forces.
        forceXyzBicycleRight(3) = 0;
        forceXyzBicycleLeft(3) = 0; 
        % Rotate crank force from bicycle reference frame to global reference
        % frame.
        forceXyzGlobalRight(1:3,iData) = rotationMatrix3d * forceXyzBicycleRight;
        forceXyzGlobalLeft(1:3,iData) = rotationMatrix3d * forceXyzBicycleLeft;

        % Set force origin as virtual marker. Convert from QTM axis to
        % OpenSim 3.3 axis. Put in column vector for .mot file format.

        % OpenSim(x) = QTM (x)
        pointXyzGlobalRight(1,iData) = virtualMarkerRight(1);
        pointXyzGlobalLeft(1,iData) = virtualMarkerLeft(1);
        % OpenSim(y) = QTM(z)
        pointXyzGlobalRight(2,iData) = virtualMarkerRight(3);
        pointXyzGlobalLeft(2,iData) = virtualMarkerLeft(3);
        % OpenSim(z) = QTM(-y)
        pointXyzGlobalRight(3,iData) = -virtualMarkerRight(2);
        pointXyzGlobalLeft(3,iData) = -virtualMarkerLeft(2);
    end
    % ----- iData -----
    
    % Write to force data file (.mot) to be used as data file in external
    % load file
    
    % Convert force to reaction force (-ve)
    reactionForceXyzGlobalRight = -forceXyzGlobalRight;
    reactionForceXyzGlobalLeft = -forceXyzGlobalLeft;
    
    % Calculate resultant force magnitude and direction just to refer to
    % later
    resultantReactionForceXyzGlobalRight = sqrt(sum(reactionForceXyzGlobalRight.^2));
    thetaReactionForceXyzGlobalRight = ...
        atan(reactionForceXyzGlobalRight(2) / reactionForceXyzGlobalRight(1));
    resultantReactionForceXyzGlobalLeft = sqrt(sum(reactionForceXyzGlobalLeft.^2));
    thetaReactionForceXyzGlobalLeft = ...
        atan(reactionForceXyzGlobalLeft(2) / reactionForceXyzGlobalLeft(1));
    
    % Convert marker coordinates from mm to SI Units (m)
    mmToM = 1000;
    pointXyzGlobalRight = pointXyzGlobalRight / mmToM;
    pointXyzGlobalLeft = pointXyzGlobalLeft / mmToM;
    
    % Create time row vector based off frame and frame rate
    nFrames = data.(trialName{1}).Frames;
    timeVector = 1/frameRate:1/frameRate:nFrames/frameRate;

    % ========== Part 2d: Visualise affect of filtering ==========
    
    % Plot signals after filtering and angle calculation
    legend1 = {'Pre-Filt','Post-Filt'};
    legend2 = {'Raw','Marker Calc.'};
    yLabel1 = 'N';
    yLabel2 = 'Deg';
    
    ax1 = subplot(3,2,1);
    plot(ax1,lTan);
    title(ax1,'ltan');
    ylabel(ax1,yLabel1)
    legend(ax1,legend1)

    ax2 = subplot(3,2,2);
    plot(ax2,rTan);
    title(ax2,'rtan');
    ylabel(ax2,yLabel1)
    legend(ax2,legend1)

    ax3 = subplot(3,2,3);
    plot(ax3,lRad);
    title(ax3,'lrad');
    ylabel(ax3,yLabel1)
    legend(ax3,legend1)

    ax4 = subplot(3,2,4);
    plot(ax4,rRad);
    title(ax4,'rrad');
    ylabel(ax4,yLabel1)
    legend(ax4,legend1)

    ax5 = subplot(3,2,[5 6]);
    plot(ax5,angleClockwiseRadians*(180/pi));
    title(ax5,'angle');
    ylabel(ax5,yLabel2)
    legend(ax5,legend2)
    pause
    savefig([folderPathSubjectResults '\' trialName{1} 'plotForceVsAngle']);
    close
      
    % ========== Part 2e: Write externalLoads.mot file ==========
    
    % Concatenate time, force and position data vertically into .mot format
    dataMot = vertcat(timeVector, reactionForceXyzGlobalRight, pointXyzGlobalRight,...
        reactionForceXyzGlobalLeft, pointXyzGlobalLeft); 
    % Change directory to subject setup folder  
    cd(folderPathSubjectSetup)
   
    % Create .mot file and open it 
    fileId = fopen([trialName{1} 'externalLoads.mot'],'w');
    % Set header names
    headers = {'time ' 'forceRightX ' 'forceRightY ' 'forceRightZ '...
    'pointRightX ' 'pointRightY ' 'pointRightZ '...
    'forceLeftX ' 'forceLeftY ' 'forceLeftZ '...
    'pointLeftX ' 'pointLeftY ' 'pointLeftZ '};
    % Write header rows into .mot file
    fprintf(fileId,'External Loads File\n');
    fprintf(fileId,'version=1\n');
    fprintf(fileId,'nRows=%d\n',size(timeVector,2));
    fprintf(fileId,'nColumns=%d\n',length(headers));
    fprintf(fileId,'Range=%d-%d seconds\n',timeVector([1 end]));
    fprintf(fileId,'endheader\n');
    str = sprintf('%s\t', headers{:});
    fprintf(fileId,'%s\t\n',str);
    % Write data into file under headers.
    fprintf(fileId,'%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n',...
        dataMot);
    fclose(fileId);
  
    % ========== Part 2f: Edit generic external load structure ==========
    
    % Edit generic external load file
    Tree.ExternalLoads.ATTRIBUTE.name = trialName{1};

    % External Loads file
    Tree.ExternalLoads.datafile = [folderPathSubjectSetup '\' trialName{1}...
        'externalLoads.mot'];

    % Motion File
    Tree.ExternalLoads.external_loads_model_kinematics_file = ...
        [folderPathSubjectResults '\' trialName{1} 'inverseKinematics.mot'] ;

    % Filter
    Tree.ExternalLoads.lowpass_cutoff_frequency_for_load_kinematics = 12;

    % ========== Part 2g: Write external load .xml file ==========

    % Set inputs for xml_write
    fileName = [folderPathSubjectSetup '\' trialName{1} 'setupExternalLoads.xml'];
    rootName = 'OpenSimDocument';
    Pref.StructItem = false;
    
    % Write .xml file
    xml_write(fileName,Tree,rootName,Pref);
    
    % ========== Part 2h: Save structure and workspace ==========

    % Save structure
    save([folderPathSubjectSetup '\' trialName{1} 'structureExternalLoads.mat'],'Tree')
    
    % Save workspace variables to use in analysis of results
    save([folderPathSubjectResults '\' trialName{1} 'workspaceExternalLoads']);
end
% ----- iFiles -----
%% End