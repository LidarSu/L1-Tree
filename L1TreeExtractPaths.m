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
% L1TreeExtractPaths.m    The function for extracting all shortest paths originating from each TLS point
%
% Version 1.0
% Latest update     27 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% fileName   Filename for point cloud data without extensions
% parameters    The parameter set
% 
% OUTPUTS:
% tlsPoints    The coordinates of TLS points and the attributes obtained based on shortest paths.
%              The first three columns record the x, y, and z coordinates of TLS points; the fourth column records the frequency of paths passing through each TLS point;
%              the fifth column records the length of the longest path passing through each TLS point; and the last column records the position of each TLS point on the longest path
% paths    The list of shorest paths. Each element represents the shortest path starting from a TLS point,consisting of the IDs of TLS points located on that path
% parameters    The updated parameter set
% ------------------------------------------------------------------------------

function [tlsPoints,paths,parameters] = L1TreeExtractPaths(fileName,parameters)

%% Load data
fileID = fopen(['./Data/',fileName,'.txt']);
tlsPoints = textscan(fileID,'%f %f %f');
x = tlsPoints{:,1};
y = tlsPoints{:,2};
z = tlsPoints{:,3};
tlsPoints = [x,y,z];
tlsPoints = unique(tlsPoints,'rows');

%% Calculate the coordinates of the base point and the base radius of the tree
tlsPoints = sortrows(tlsPoints,3);
bottomPoints = tlsPoints(1:fix(size(tlsPoints,1)/100),:);
minX = quantile(bottomPoints(:,1),0.05);
maxX = quantile(bottomPoints(:,1),0.95);
minY = quantile(bottomPoints(:,2),0.05);
maxY = quantile(bottomPoints(:,2),0.95);
minZ = quantile(bottomPoints(:,3),0.05);
maxZ = quantile(bottomPoints(:,3),0.95);
selectIDs = find((bottomPoints(:,1) >= minX) & (bottomPoints(:,1) <= maxX) & ...
    (bottomPoints(:,2) >= minY) & (bottomPoints(:,2) <= maxY) & ...
    (bottomPoints(:,3) >= minZ) & (bottomPoints(:,3) <= maxZ));
bottomPoints = bottomPoints(selectIDs,:);
oriParameters = [0,0,0];
F = @(k) (bottomPoints(:,1)-k(1)).^2 + (bottomPoints(:,2) - k(2)).^2 - k(3)^2;    % Calculate coordinates using the Gauss-Newton method
optParameters = lsqnonlin(F,oriParameters);
tlsPoints = [tlsPoints;[optParameters(1),optParameters(2),minZ]];
parameters.basePoint = [optParameters(1),optParameters(2),minZ];
parameters.baseRadius = abs(optParameters(3));
parameters.initialSearchRange = parameters.baseRadius*0.1;
parameters.searchRange = parameters.baseRadius*0.1;

%% Calculate shortest paths
[kDistances,kIDs] = pdist2(tlsPoints,tlsPoints,'euc','Smallest',11);
kDistances = kDistances(2:end,:);
kIDs = kIDs(2:end,:);
startPointIDs = repmat(linspace(1,size(kDistances,2),size(kDistances,2)),size(kDistances,1),1);
startPointIDs = reshape(startPointIDs,1,size(kDistances,1)*size(kDistances,2));
targetPointIDs = reshape(kIDs,1,size(kDistances,1)*size(kDistances,2));
weightMatrix = reshape(kDistances,1,size(kDistances,1)*size(kDistances,2));
treeGraph = graph(startPointIDs,targetPointIDs,weightMatrix);
paths = shortestpathtree(treeGraph,'all',size(tlsPoints,1),'OutputForm','cell');
paths = paths(1:(end-1),:);
attrMatrix = zeros(size(tlsPoints,1),3);
for i = 1:1:size(paths,1)
    tempIDs = paths{i};
    tempLocations = linspace(1,length(tempIDs),length(tempIDs));
    attrMatrix(tempIDs,1) = attrMatrix(tempIDs,1) + 1;
    selectIDs = find(attrMatrix(tempIDs,2) < length(tempIDs));
    if ~isempty(selectIDs)
        attrMatrix(tempIDs(selectIDs),2) = length(tempIDs);
        attrMatrix(tempIDs(selectIDs),3) = tempLocations(selectIDs);
    end
end
tlsPoints = [tlsPoints,attrMatrix];

end