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
% L1TreeInitialization.m    The initialization function of the L1-Tree algorithm
% 
% Version 1.0
% Latest update     27 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% OUTPUTS:
% parameters.leafOn    This parameter indicates whether the point cloud contains leaf points
% parameters.woodProp    This parameter is related the number of points classified as wood. The larger the parameter, the more points will be classified as wood
% parameters.basePoint    The coordinates of the circle center at the base of the tree
% parameters.baseRadius    The base radius of the tree
% parameters.initialSearchRange    Initial searching range, typically set as a multiple of the base radius
% parameters.minMu    Penalty coefficient
% parameters.maxMu    Penalty coefficient
% parameters.maxIterationNumFixedRange    Number of iterations for calculating skeleton point coordinates with a constant searching range
% parameters.iterationSpeed    Incremental multiplier for the searching range
% parameters.nearThreshold    Minimum spacing between skeleton points
% parameters.sigmaThreshold    Directional degree (sigma) threshold required for identifying a skeleton point as the initial point of a skeleton segment
% parameters.kSmoothSigma    Number of neighboring skeleton points used for smoothing sigma values
% parameters.sigmaTracingThreshold    Directional degree threshold necessary for the iterative search of the next point in a skeleton segment
% parameters.distanceTracingThreshold    Distance threshold necessary for the iterative search of the next point in a skeleton segment
% parameters.angleTracingThreshold    Angle threshold necessary for the iterative search of the next point in a skeleton segment
% parameters.kTracingThreshold    Maximum number of candidates for the next point of a skeleton segment
% parameters.branchSizeThreshold    Length threshold for a skeleton segment or a branch, counted by the minimum number of skeleton points required to form a segment or branch
% parameters.minDistanceThreshold    Distance threshold for creating connections without constraints
% parameters.maxDistanceThreshold    Maximum length of a connection
% parameters.minAngleThreshold    Minimum angle threshold for examining the validity of a connection
% parameters.maxAngleThreshold    Maximum angle threshold for examining the validity of a connection
% parameters.densityFactor    Coefficient of point density used to examine the validity of a connection
% parameters.resolution    Spacing between skeleton points added to a connection
% parameters.windowSize    Sliding window size used for calculating linear density
% parameters.moveStep    Moving step of the sliding window used for calculating linear density
% parameters.smoothAngleThreshold    Angle threshold for smoothing a branch
% ------------------------------------------------------------------------------

function parameters = L1TreeInitialization()

parameters.leafOn = false;
parameters.woodProp = 0.03;
parameters.basePoint = 0;
parameters.baseRadius = 0;
parameters.initialSearchRange = 0;
parameters.searchRange = 0;
parameters.minMu = 0.15;
parameters.maxMu = 0.35;
parameters.maxIterationNumFixedRange = 20;
parameters.iterationSpeed = 0.4;
parameters.nearThreshold = 0.01;
parameters.sigmaThreshold = 0.9;
parameters.kSmoothSigma = 6;
parameters.sigmaTracingThreshold = 0.8;
parameters.distanceTracingThreshold = 0.4;
parameters.angleTracingThreshold = 25;
parameters.kTracingThreshold = 12;
parameters.branchSizeThreshold = 3;

parameters.minDistanceThreshold = 0.1;
parameters.maxDistanceThreshold = 0.8;
parameters.minAngleThreshold = 10;
parameters.maxAngleThreshold = 120;
parameters.densityFactor = 1/3;
parameters.resolution = 0.03;
parameters.windowSize = 0.1;
parameters.moveStep = 0.01;
parameters.smoothAngleThreshold = 45;

end