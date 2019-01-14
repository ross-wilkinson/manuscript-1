% ########## Setup & Run Point Kinematics ##########
% 
% Description
% ---------------
% This script writes a subject specific setup file to analyze point
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

% Flag
Tree.AnalyzeTool.AnalysisSet.objects.PointKinematics.on = 'true';

% Set joint specific variables before looping through
jointList = {'hip','knee','ankle'};
nJoints = numel(jointList);

%% ========== Part 3: Loop through joints ==========
% Progress counter
j = 0;

% ----- iJoints -----
for iJoints = 1:nJoints
    if strcmp(jointList{iJoints},'hip')
        % child
        bodyName = 'femur_r'; 
        % reference system
        relativeBody = 'ground';
        pointName = 'CorHipRight';
    elseif strcmp(jointList{iJoints},'knee')
        bodyName = 'tibia_r'; 
        relativeBody = 'ground';
        pointName = 'CorKneeRight';
    else
        bodyName = 'talus_r';
        relativeBody = 'ground';
        pointName = 'CorAnkleRight';
    end
    
    % Body to analyze
    Tree.AnalyzeTool.AnalysisSet.objects.PointKinematics.body_name = bodyName; 
    
    % Ref body
    Tree.AnalyzeTool.AnalysisSet.objects.PointKinematics. ...
        relative_to_body_name = relativeBody;
    
    % Point name
    Tree.AnalyzeTool.AnalysisSet.objects.PointKinematics.point_name = pointName;
    
    % Location within child body to track
    Tree.AnalyzeTool.AnalysisSet.objects.PointKinematics.point = [0 0 0]; 

    %% ========== Part 3a: Loop through each subject file ==========
    
    % Display progress of script
    h = waitbar(0,...
        ['Running ' subjectName ' Point Kinematics.']);
    
    % ----- iFiles -----
    for iFiles = 1:nFiles
                
        % Update progress
        j = j+1;
        progress = j/(nFiles*nJoints);
        j = j+1;
        n = num2str(nFiles*nJoints);

        % Compute percentage of completion
        waitbar(progress,h,...
            ['Running ' subjectName ' Point Kinematics. File ' j ' of ' n])
    
        cd(folderPatghSubjectData)
        data = load(subjectFileList(iFiles).name);
        trialName = fieldnames(data); 
        
        % Get times
        frameRate = data.(trialName{1}).FrameRate;
        nFrames = data.(trialName{1}).Frames;
        startFrame = data.(trialName{1}).StartFrame;
        endFrame = nFrames + startFrame - 1;
        startTime = startFrame / frameRate;
        endTime = endFrame / frameRate;
        
        % ========== Part 3b: File edits to generic Analyze structure ==========
        
        % Set initial time
        Tree.AnalyzeTool.initial_time = startTime; 
        
        % Set final time
        Tree.AnalyzeTool.final_time = endTime;
        
        % Analysis name
        Tree.AnalyzeTool.ATTRIBUTE.name = trialName{1};
        
        % Start time
        Tree.AnalyzeTool.AnalysisSet.objects.PointKinematics.start_time = ...
            startTime;
        
        % End time
        Tree.AnalyzeTool.AnalysisSet.objects.PointKinematics.end_time = endTime; 
        
        % Coordinates file
        Tree.AnalyzeTool.coordinates_file = ...
            [folderPathSubjectResults '\' trialName{1} 'inverseKinematics.mot']; 
        
        % ========== Part 3c: Write .xml setup file ==========
        
        % Set inputs for xml_write  
        fileName = [folderPathSubjectSetup '\' ...
            trialName{1} 'setupAnalyzePointKinematics' pointName '.xml'];
        rootName = 'OpenSimDocument';
        Pref.StructItem = false;
        
        % Write .xml file
        xml_write(fileName,Tree,rootName,Pref); %write .xml file

        % ========== Part 3d: Save structure ==========
        
        % Save structure
        save([folderPathSubjectSetup '\' ...
            trialName{1} 'structureAnalyzePointKinematics' pointName '.mat'],...
            'Tree')
        
        % ========== Part 3e: Run Analyze Tool ==========
        
        % Include double inverted commas so that the command line gets rid of
        % spaces in path name.
        analyzeFile = ['"' fileName '"'];
        command = ['analyze -S ' analyzeFile]; 
        system(command);
    end
    % ----- iFiles -----
end
% ----- iJoints -----

%% End