% ########## Setup & Run Inverse Dynamics ##########
% 
% Description
% ---------------
% This script writes a subject specific setup file to be used with the
% Inverse Dynamics Tool in OpenSim 3.3 and then runs the ID Tool through the
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
% of Queensland) MatlabOpensimPipelineTools 'setup_InverseDynamics.m'
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

% Load structure of generic ID file
load('D:\opensim\structureInverseDynamics.mat')

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

% Find subject's model file
cd(folderPathSubjectResults)
modelFileList = dir('*.osim');
modelFile = [modelFileList.folder '\' modelFileList.name];

% InverseDynamicsTool -> Directories
Tree.InverseDynamicsTool.results_directory = folderPathSubjectResults;

% InverseDynamicsTool -> Model file
Tree.InverseDynamicsTool.model_file = modelFile;

%% ========== Part 2: Loop through each subject file ==========

% Display progress of script
h = waitbar(0,...
    ['Running ' subjectName ' Inverse Dynamics.']);

% ----- iFiles -----
for iFiles = 1:nFiles
    
    % Update progress
    progress = iFiles/nFiles;
    j = num2str(iFiles);
    n = num2str(nFiles);
    
    % Compute percentage of completion
    waitbar(progress,h,...
        ['Running ' subjectName ' Inverse Dynamics. File ' j ' of ' n])
    
    cd(folderPathSubjectData)
    data = load(subjectFileList(iFiles).name);
    trialName = fieldnames(data); 
    
    % ========== Part 2a: Edit generic ID structure ==========
    
    % InverseDynamicsTool
    Tree.InverseDynamicsTool.ATTRIBUTE.name = trialName{1};
    
    % Get times
    frameRate = data.(trialName{1}).FrameRate;
    nFrames = data.(trialName{1}).Frames;
    startFrame = data.(trialName{1}).StartFrame;
    endFrame = nFrames + startFrame - 1;
    startTime = startFrame / frameRate;
    endTime = endFrame / frameRate;
    
    % InverseDynamicsTool -> Time range
    Tree.InverseDynamicsTool.time_range = [startTime endTime];
    
    % InverseDynamicsTool -> External Loads File
    Tree.InverseDynamicsTool.external_loads_file = ...
        [folderPathSubjectSetup '\' trialName{1} 'setupExternalLoads.xml'];

    % InverseDynamicsTool -> Motion File
    Tree.InverseDynamicsTool.coordinates_file = ...
        [folderPathSubjectResults '\' trialName{1} 'inverseKinematics.mot'];

    % InverseDynamicsTool -> Output File
    Tree.InverseDynamicsTool.output_gen_force_file = [trialName{1} 'inverseDynamics.sto'];

    % ========== Part 2b: Write .xml setup file ==========
    
    % Set inputs for xml_write
    fileName = [folderPathSubjectSetup '\' trialName{1} 'setupInverseDynamics.xml'];
    rootName = 'OpenSimDocument';
    Pref.StructItem = false;
    
    % Write .xml file
    xml_write(fileName,Tree,rootName,Pref);

    % ========== Part 2c: Save structure ==========
    
    % Save structure
    save([folderPathSubjectSetup '\' trialName{1} 'structureInverseDynamics.mat'],'Tree')
    
    % ========== Part 2d: Run ID Tool ==========
    
    % Include double inverted commas so that the command line gets rid of
    % spaces in path name.
    inverseDynamicsFile = ['"' fileName '"'];
    command = ['id -S ' inverseDynamicsFile];
    system(command);
end
% ----- iFiles -----

%% End