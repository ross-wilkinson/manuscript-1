% ########## Setup & Run Inverse Kinematics ##########
% 
% Description
% ---------------
% This script writes a subject specific setup file to be used with the
% Inverse Kinematics Tool in OpenSim 3.3 and then runs the IK Tool through the
% command line.
%
% Pre-requisites
% --------------- 
% To execute this script you will need to:
%   1. Download the MATLAB-OPENSIM INTERFACE TOOLBOX which can be accessed
%   here: https://simtk.org/projects/matlab_tools
%   2. Create a generic structure file which has many of the structure
%   fields pre-filled with default settings. This can be done by running
%   the 'writeSetupInverseDynamicsXml.m' script, which is available within
%   my other GitHub repo: https://github.com/uqrwilk/opensimModelling
%
% Acknowledgement
% ----------------
% This script is merely an adaptation of Glen Lichtwark's (The University
% of Queensland) MatlabOpensimPipelineTools 'setup_InverseKinematics.m'
% function. The functions available within the above mentioned toolbox are
% far more efficient and robust than my scripts, however I made these
% scripts as a learning process for myself.
% 
% Written by Ross Wilkinson, The University of Queensland,
% (ross.wilkinson_at_uqconnect.edu.au)

%% Initialise
clear;clc;close all
import org.opensim.modeling.*

%% ========== Part 1: Get and set subject data and folders ==========

% Load structure of generic IK file
load('D:\opensim\structureInverseKinematics.mat')

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

% Find scaled subject model file
cd(folderPathSubjectResults)
modelFileList = dir('*.osim');
modelFile = [modelFileList.folder '\\' modelFileList.name];

%% ========== Part 2: Loop through each subject file ==========

% Display progress of script
h = waitbar(0,...
    ['Running ' subjectName ' Inverse Kinematics.']);

% ----- iFiles -----
for iFiles = 1:nFiles
    
    % Update progress
    progress = iFiles/nFiles;
    j = num2str(iFiles);
    n = num2str(nFiles);
    
    % Fill percentage of waitbar
    waitbar(progress,h,...
        ['Running ' subjectName ' Inverse Kinematics. File ' j ' of ' n])
    
    cd(folderPathSubjectData)
    data = load(subjectFileList(iFiles).name);
    trialName = fieldnames(data);
    
    % ========== Part 2a: Rotate coordinate system ==========
    
    % Rotate markers from global to bicycle coordinate system
    markerData = data.(trialName{1}).Trajectories.Labeled.Data(:,1:3,:);
    markerLabels = categorical(data.(trialName{1}).Trajectories.Labeled.Labels);

    % Re-shape marker_data for writing .trc file. One marker per page.
    markerData = permute(markerData,[2 3 1]);

    % Assign bike markers to variables
    B2 = markerData(:,1,markerLabels == 'B2');
    B3 = markerData(:,1,markerLabels == 'B3');

    % Use bike markers to create bicycle coordinate system
    dist = B2' - B3'; % transpose for vertcat
    
    % New y-axis in bicycle reference frame. Divide by Euclidean length of the
    %vector.
    yAxisBicycle = dist / norm(dist);
    
    % Use global z-axis to calculate x-axis in bicycle reference frame.
    zAxisBicycle = [0 0 1];
    
    % Take cross product of y(bicycle) and z(global) to find new x-axis in
    % bicycle reference frame.
    xAxisBicycle = cross(yAxisBicycle,zAxisBicycle);
    
    % Now use x and y-axis in bicycle reference frame to calculate new z-axis.
    zAxisBicycle = cross(xAxisBicycle,yAxisBicycle);
    
    % Set bicycle unit vectors in rotation matrix
    rotationMatrix3d = [xAxisBicycle;yAxisBicycle;zAxisBicycle];

    % Rotate crank force from bicycle reference frame to global reference
    % frame.
    markerDataBicycle = zeros(size(markerData));
    
    % ----- i -----
    for i = 1:size(markerData,3)
        markerDataBicycle(:,:,i) = rotationMatrix3d * markerData(:,:,i);
    end
    % ----- i -----

    % Change axes to match OpenSim 3.3
    markersToOpenSim = zeros(size(markerDataBicycle));
    % x = x
    markersToOpenSim(1,:,:) = markerDataBicycle(1,:,:);
    % y = z
    markersToOpenSim(2,:,:) = markerDataBicycle(3,:,:);
    % z = -y
    markersToOpenSim(3,:,:) = -markerDataBicycle(2,:,:);

    % ========== Part 2b: Write .trc file ==========

    % Add .trc file ext to filename
    fileNameTrc = [trialName{1} '.trc'];

    % Open the file
    fileId = fopen([folderPathSubjectData '\' fileNameTrc],'w');

    fprintf(fileId,'PathFileType\t4\t(X/Y/Z)\t %s\n',fileNameTrc);
    fprintf(fileId,...
        'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
    fprintf(fileId,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n',...
        data.(trialName{1}).FrameRate,data.(trialName{1}).FrameRate,...
        data.(trialName{1}).Frames,...
        data.(trialName{1}).Trajectories.Labeled.Count,'m',...
        data.(trialName{1}).FrameRate,data.(trialName{1}).StartFrame,...
        data.(trialName{1}).Frames);
    % create header4 using marker labels
    header4 = 'Frame#\tTime\t';
    % create header5 using XYZ and marker numbers
    header5 = '\t\t';
    formatText = '%i\t%2.4f\t';
    startFrame = data.(trialName{1}).StartFrame;
    endFrame = data.(trialName{1}).Frames + startFrame - 1;
    nFrames = startFrame:endFrame;
    frameRate = data.(trialName{1}).FrameRate;
    time = nFrames/frameRate;
    dataOutput = [nFrames; time];
    
    % ----- mm -----
    for mm = 1:length(markerLabels)
        header4 = [header4 char(markerLabels(mm)) '\t\t\t'];
        header5 = [header5 'X' num2str(mm) '\t' 'Y' num2str(mm) '\t'...
            'Z' num2str(mm) '\t'];
        formatText = [formatText '%f\t%f\t%f\t'];
        dataOutput = [dataOutput;markersToOpenSim(:,:,mm)/1000];
    end
    % ----- mm -----
    
    header4 = [header4 '\n'];
    header5 = [header5 '\n'];
    formatText = [formatText '\n'];

    fprintf(fileId,header4);
    fprintf(fileId,header5);
    fprintf(fileId,formatText,dataOutput);
    fclose(fileId);
    
    % ========== Part 2c: Edit IK setup file ==========
    
    % Assign .trc file to Marker File
    markerFile = ([folderPathSubjectData '\' fileNameTrc]);
    % Edit Inverse Kinetmatics Tool
    Tree.InverseKinematicsTool.ATTRIBUTE.name = trialName{1};
    % Edit Inverse Kinetmatics Tool -> Directories
    Tree.InverseKinematicsTool.results_directory = folderPathSubjectResults;
    % Edit Inverse Kinetmatics Tool -> Model file
    Tree.InverseKinematicsTool.model_file = modelFile;
    % Edit Inverse Kinetmatics Tool -> Marker File
    Tree.InverseKinematicsTool.marker_file = markerFile;
    % Edit Inverse Kinematics Tool -> Time Range
    startTime = startFrame/data.(trialName{1}).FrameRate;
    endTime = endFrame/data.(trialName{1}).FrameRate;
    Tree.InverseKinematicsTool.time_range = [startTime endTime];
    % Edit Inverse Kinetmatics Tool -> Output file
    Tree.InverseKinematicsTool.output_motion_file = ...
        [folderPathSubjectResults '\' trialName{1} 'inverseKinematics.mot'];
    
    % ========== Part 2d: Setup and write .xml file ==========
    
    % Set inputs for xml_write
    fileName = [folderPathSubjectSetup '\' trialName{1}...
        'setupInverseKinematics.xml'];
    rootName = 'OpenSimDocument';
    Pref.StructItem = false;
    
    % Write .xml file
    xml_write(fileName,Tree,rootName,Pref);

    % ========== Part 2e: Save structure ==========
    save([folderPathSubjectSetup '\' trialName{1}...
        'structureInverseKinematics.mat'],'Tree')

    % ========== Part 2f: Run IK Tool ==========
    
    % Include double inverted commas so that the command line gets rid of
    % spaces in path name.
    inverseKinematicsFile = ['"' fileName '"'];
    command = ['ik -S ' inverseKinematicsFile];
    system(command);
end
% ----- iFiles -----

%% End