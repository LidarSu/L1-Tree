% This File is part of the L1-Tree algorithm
% 
% L1-Tree is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% 
% If you use this file, please cite the article entitled "L1-Tree: 
% A novel algorithm for constructing 3D tree models and estimating branch architectural 
% traits using terrestrial laser scanning data"

% ------------------------------------------------------------------------------
% Main.m    The main function of the L1-Tree algorithm
%
% Version 1.0
% Latest update     26 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% fileName    Filename for point cloud data without extensions
% parameter.leafOn    Whether the input point cloud contains leaf points
% parameters.skeletonProp    Activate this parameter when parameter.leafOn is TRUE. This parameter regulates the bulit-in leaf-wood
%                            seperation algorithm. The larger the value of this parameter, the more points will be classified as wood points
% 
% OUTPUTS:
% fileName_L1Median_SeedPoints.txt    The output file of the bulit-in L1-Median algorithm, which records seed point coordinates
% fileName_L1Median_Branches.txt    The output file of the bulit-in L1-Median algorithm, which records the iniitial skeleton
% fileName_L1Tree_3DModel.txt    The output file of the L1-Tree algorithm, which can be used as an input file for drawing the 3D tree model in Python
% fileName_L1Tree_Centrifugal.mat    The output MATLAB file of the L1-Tree algorithm, which records the attributes of each branch  (e.g., ID, code, length and radius) under the centrifugal coding strategy
% fileName_L1Tree_Gaaliche.mat    The output MATLAB file of the L1-Tree algorithm, which records the attributes of each branch  (e.g., ID, code, length and radius) under Gaaliche's coding strategy
% ------------------------------------------------------------------------------

%% Initialize
clc,clear;
% Load the Python environment
% pyversion('C:\Users\erbao\anaconda3\envs\Draw3d\python.exe');    % Please eusure that the Python is version 3.7 and contains required modules, including pandas, matplotlib, PIL, scipy, and mayavi
% pyModule = py.importlib.import_module('Draw');
% py.importlib.reload(pyModule);
% Load default parameter values
parameters = L1TreeInitialization();    % You can customise the parameter values by modifying the L1TreeInitialization.m file
% parameters.leafOn = true;    % If parameters.leafOn = true, the built-in leaf-wood seperation algorithm will be activated and a suitable parameters.woodProp value needs to be provided
% parameters.woodProp = 0.03;
% Input file name
fileName = 'Tree_1';    % Please ensure the Tree_1.txt file is in the Data folder

%% Select seed points
[tlsPoints,paths,parameters] = L1TreeExtractPaths(fileName,parameters);
[seedPoints,tlsPoints] = L1TreeExtractSeedPoints(tlsPoints,paths,parameters);

%% Run the built-in L1-Median algorithm and extract the initial skeleton
globalBar = waitbar(0,'L_1 Median: Start');    % The globalBar is used to show the step being execuated by the L1-Tree algorithm
iterationNum = 0;
skeletonBranches = [];
skeletonPoints = seedPoints;
maxIterationNum = ceil(log(parameters.baseRadius*2/parameters.searchRange)/log(parameters.iterationSpeed + 1))*parameters.maxIterationNumFixedRange;
tlsPoints = L1MedianCalculateDensity(tlsPoints,parameters,globalBar,iterationNum,maxIterationNum);
while parameters.searchRange <= parameters.baseRadius*2
    iterationNumFixedRange = 0;
    averageMovingDistance = 0;
    while iterationNumFixedRange < parameters.maxIterationNumFixedRange
        localBar = waitbar(0,'One iteration: start','Position',[100,100,270,60]);    % The localBar is used to show the information about the current iteration, including the search range, the iteration number, and the average moving distance
        skeletonPoints = L1MedianSkeletonNeighbors(skeletonPoints,tlsPoints,parameters,localBar);
        skeletonPoints = L1MedianCalculateSigma(skeletonPoints,parameters,localBar);
        skeletonPoints = L1MedianCalculateAverageTerm(skeletonPoints,parameters,localBar);
        skeletonPoints = L1MedianCalculateRepulsionTerm(skeletonPoints,parameters,localBar);
        [skeletonPoints,averageMovingDistance] = L1MedianUpdateSkeletonPoints(skeletonPoints,parameters,localBar);
        iterationNumFixedRange = iterationNumFixedRange + 1;
        iterationNum = iterationNum + 1;
        close(localBar);
        showInfo = ['L_1 Median: ',num2str(parameters.searchRange),' ',num2str(iterationNumFixedRange),' ',num2str(averageMovingDistance)];
        waitbar(iterationNum/maxIterationNum,globalBar,showInfo);
    end
    skeletonPoints = L1MedianRemoveSkeletonPoints(skeletonPoints,parameters,globalBar,iterationNum,maxIterationNum);
    skeletonPoints = L1MedianSmoothSigma(skeletonPoints,parameters,globalBar,iterationNum,maxIterationNum);
    [newBranches,skeletonPoints] = L1MedianSearchNewBranches(skeletonPoints,parameters,globalBar,iterationNum,maxIterationNum);
    skeletonBranches = [skeletonBranches,newBranches];
    parameters.searchRange = parameters.searchRange*(1 + parameters.iterationSpeed);
end
close(globalBar);

%% Optimize the skeleton and construct 3D tree model
globalBar = waitbar(0,'MST optimize: Start');
[treeModel,tlsPoints] = L1TreeConnect(skeletonBranches,skeletonPoints,tlsPoints,parameters,globalBar);
[treeModel,tlsPoints] = L1TreePruneBranches(treeModel,tlsPoints,parameters,globalBar);
[treeModel,tlsPoints] = L1TreeAdjustBifurcation(treeModel,tlsPoints,parameters,globalBar);
[treeModel,tlsPoints] = L1TreeSmoothBranches(treeModel,tlsPoints,parameters,globalBar);
[treeModel,tlsPoints] = L1TreeExtendSkeleton(treeModel,tlsPoints,parameters,globalBar);
[treeModelCentrifugal,treeModelGaaliche] = L1TreeConstructModel(fileName,treeModel,tlsPoints,globalBar);
close(globalBar);

%% Draw 3D tree model
% textureFileName = './Data/bark.jpg';
drawFileName = L1TreeExportModel(fileName,treeModelCentrifugal);
% py.Draw.Draw3DTreeModel(drawFileName,textureFileName);