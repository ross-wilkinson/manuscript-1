% ########## Setup Scale ##########
% 
% Description
% -----------
% This script writes a subject specific setup file to be used with the 
% Scale Tool in OpenSim 3.3.
%
% Pre-requisites
% -------------- 
% 1. To execute this script you will need to download the MATLAB-OPENSIM
% INTERFACE TOOLBOX which can be accessed here:
% https://simtk.org/projects/matlab_tools
% 
% 2. Create a generic structure file which has many of the structure
%   fields pre-filled with default settings. This can be done by running
%   the 'writeSetupInverseDynamicsXml.m' script, which is available within
%   my other GitHub repo: https://github.com/uqrwilk/opensimModelling
%
% Acknowledgement
% ----------------
% This script is merely an adaptation of Glen Lichtwark's (The University
% of Queensland) MatlabOpensimPipelineTools 'setup_scale.m'
% function. The functions available within the above mentioned toolbox are
% far more efficient and robust than my scripts, however I made these
% scripts as a learning process for myself.
% 
% Written by Ross Wilkinson, The University of Queensland,
% (ross.wilkinson_at_uqconnect.edu.au)

%% Initialise
clear;clc;close all

%% ========== Part 1: Get and set subject data and folders ==========

% Load structure of generic scale file
load('D:\opensim\structureScale.mat')

% Input subject number
subjectNo = input('Subject Number? ','s');

% Create subject name
subjectName = ['subject' subjectNo];

% Set filepath for subject setup and results files
folderPathExp = 'D:\exp01';
folderPathSubSetup = [folderPathExp '\' subjectName '\setup'];
folderPathSubResults = [folderPathExp '\' subjectName '\results'];
folderPathSubData = [folderPathExp '\' subjectName '\data'];

% Find subject details in data table
dataTable = readtable([folderPathExp '\dataTable.xlsx']);
iSubject = find(contains(dataTable.subject,subjectName));
subjectMass = dataTable.mass(iSubject);
subjectHeight = dataTable.height(iSubject);
testDate = datenum(dataTable.testDate(iSubject));
dob = datenum(dataTable.dob(iSubject));
subjectAge = str2double(datestr(testDate-dob,11));

%% ========== Part 2: Convert .c3d file to .trc ==========

fileExtC3d = '_static.c3d';
fileNameC3d = [folderPathSubData '\\' subjectName fileExtC3d]; 
data = btk_c3d2trc(fileNameC3d);
fileNameTrc = regexprep(fileNameC3d,'.c3d','.trc');

%% ========== Part 3: Create structure ==========

% Edit ScaleTool structure
Tree.ScaleTool.ATTRIBUTE.name = subjectName;
Tree.ScaleTool.mass = subjectMass;
Tree.ScaleTool.height = subjectHeight;
Tree.ScaleTool.age = subjectAge;

% Edit ScaleTool -> GenericModelMaker
Tree.ScaleTool.GenericModelMaker.ATTRIBUTE.name = subjectName;

% Edit ScaleTool -> ModelScaler
Tree.ScaleTool.ModelScaler.ATTRIBUTE.name = subjectName;
Tree.ScaleTool.ModelScaler.marker_file = fileNameTrc;

% ScaleTool -> MarkerPlacer
Tree.ScaleTool.MarkerPlacer.ATTRIBUTE.name = subjectName;
Tree.ScaleTool.MarkerPlacer.marker_file = fileNameTrc;
Tree.ScaleTool.MarkerPlacer.output_model_file = ...
    [folderPathSubResults '\' subjectName 'modelScaled.osim'];
Tree.ScaleTool.MarkerPlacer.output_motion_file = ...
    [folderPathSubResults '\' subjectName 'static.mot'];
Tree.ScaleTool.MarkerPlacer.output_marker_file = ...
[folderPathSubResults '\' subjectName 'markersScaled.xml'];

%% ========== Part 4: Convert structure to .xml file ==========

% Set inputs for xml_write
fileNameScaleXml = [folderPathSubSetup '\' subjectName 'setupScale.xml'];
rootName = 'OpenSimDocument';
Pref.StructItem = false;

% Write .xml file
xml_write(fileNameScaleXml,Tree,rootName,Pref);

%% ========== Part 5: Save structure ==========

save([folderPathSubSetup '\' subjectName 'structureScale.mat'],'Tree')
disp('Done.')

%% END