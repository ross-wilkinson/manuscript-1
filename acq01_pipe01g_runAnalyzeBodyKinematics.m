% ########## Setup & Run Body Kinematics ##########
% 
% Description
% ---------------
% This script writes a subject specific setup file to analyze body
% kinematics within the Analyze Tool in OpenSim 3.3 and then runs the
% Analyze Tool through the command line.
%
% Pre-requisites
% --------------- 
% To execute this script you will need to:
%   1. Download the MATLAB-OPENSIM INTERFACE TOOLBOX which can be accessed
%   here: https://simtk.org/projects/matlab_tools
%   2. Create a generic structure file which has many of the structure
%   fields pre-filled with default settings. This can be done by running
%   the 'writeSetupAnalyzeXml.m' script, which is available within
%   my other GitHub repo: https://github.com/uqrwilk/opensimModelling
%
% Acknowledgement
% ----------------
% This script is merely an adaptation of Glen Lichtwark's (The University
% of Queensland) MatlabOpensimPipelineTools 'setup_MuscleAnalysis.m'
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

% Load structure of generic Analyze file
load('D:\opensim\structureAnalyze.mat')

% Input subject number
subjectNo = input('Subject Number? ','s'); 

% Create subject name
subjectName = ['subject' subjectNo]; 

% Set filepath for subject setup and results files
folderPathExperiment = 'D:\exp01'; 
folderPathSubjectSetup = [folderPathExperiment '\' subjectName '\setup'];
folderPathSubjectResults = [folderPathExperiment '\' subjectName '\results'];
folderPatghSubjectData = [folderPathExperiment '\' subjectName '\data'];

% Find subject's exported .mat files
cd(folderPatghSubjectData)
subjectFileList = dir('*.mat'); 
nFiles = size(subjectFileList,1);

% Find subject's model file
cd(folderPathSubjectResults)
modelFileList = dir('*.osim'); 
modelFile = [modelFileList.folder '\' modelFileList.name];

%% ========== Part 2: Initial edit of generic structure ==========

% Model file
Tree.AnalyzeTool.model_file = modelFile; 

% Results directory
Tree.AnalyzeTool.results_directory = folderPathSubjectResults; 

% turn on Body Kinematics
Tree.AnalyzeTool.AnalysisSet.objects.BodyKinematics.on = 'true';

% Set coordinates to report
Tree.AnalyzeTool.AnalysisSet.objects.BodyKinematics.bodies = 'center_of_mass';

%% ========== Part 3: Loop through each subject file ==========

% Display progress of script
h = waitbar(0,...
    ['Running ' subjectName ' Body Kinematics.']);

% ----- iFiles -----
for iFiles = 1:nFiles
    
    % Update progress
    progress = iFiles/nFiles;
    j = num2str(iFiles);
    n = num2str(nFiles);
    
    % Compute percentage of completion
    waitbar(progress,h,...
        ['Running ' subjectName ' Body Kinematics. File ' j ' of ' n])
    
    cd(folderPathSubjectData)
    data = load(subjectFileList(iFiles).name);
    trialName = fieldnames(data);
    
    % Get times
    frameRate = data.(trialName{1}).FrameRate;
    nFrames = data.(trialName{1}).Frames;
    startFrame = data.(trialName{1}).StartFrame;
    endFrame = nFrames + startFrame - 1;
    startTime = startFrame / frameRate;
    endTime = endFrame / frameRate;
    
    % ========== Part 3a: Edit generic Analyze structure ==========
    
    % Intial time
    Tree.AnalyzeTool.initial_time = startTime;
    
    % End time
    Tree.AnalyzeTool.final_time = endTime;
    
    % Analysis name
    Tree.AnalyzeTool.ATTRIBUTE.name = trialName{1};
    
    % Start time
    Tree.AnalyzeTool.AnalysisSet.objects.BodyKinematics.start_time = startTime; 
    
    % End time
    Tree.AnalyzeTool.AnalysisSet.objects.BodyKinematics.end_time = endTime; 
    
    % Coordinates file
    Tree.AnalyzeTool.coordinates_file = ...
        [folderPathSubjectResults '\' trialName{1} 'inverseKinematics.mot']; 
    
    % ========== Part 3b: Write .xml setup file ==========
    
    % Set inputs for xml_write
    fileName = [folderPathSubjectSetup '\' trialName{1} 'setupAnalyzeBodyKinematics.xml'];
    rootName = 'OpenSimDocument';
    Pref.StructItem = false;
    
    % Write .xml file
    xml_write(fileName,Tree,rootName,Pref);
    
    % ========== Part 3c: Save structure ==========
    
    % Save structure
    save([folderPathSubjectSetup '\' trialName{1} 'structureAnalyzeBodyKinematics.mat'],'Tree');

    % ========== Part 3d: Run Analyze Tool ==========
    
    % Include double inverted commas so that the command line gets rid of
    % spaces in path name.
    analyzeFile = ['"' fileName '"'];
    command = ['analyze -S ' analyzeFile];
    system(command);
end
% ----- iFiles -----

%% End