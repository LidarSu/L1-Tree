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
% L1TreeExtractSeedPoints.m    The function for extract seed points
%
% Version 1.0
% Latest update     27 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% tlsPoints    The coordinates of TLS points and the attributes obtained based on shortest paths.
%              The first three columns record the x, y, and z coordinates of TLS points; the fourth column records the frequency of paths passing through each TLS point;
%              the fifth column records the length of the longest path passing through each TLS point; and the last column records the position of each TLS point on the longest path
% paths    The list of shorest paths. Each element represents the shortest path starting from a TLS point,consisting of the IDs of TLS points located on that path
% parameters    The updated parameter set
% 
% OUTPUTS:
% seedPoints    The struct that records the attributes of skeleton points.
%               The attributes include the ID of a skeleton point (ID); the coordinates of a skeleton point (P);
%               the length of the longest path passing through a skeleton point (PathLength); the location of a skeleton point on the longest path passing through it (Location);
%               the flag indicating whether a skeleton point is fixed or not (IsFix); the flag indicating whether a skeleton point needs to be ignored or not (IsIgnore);
%               the flag indicating whether a skeleton point has become a branch point or not (IsBranch); the IDs of neighboring skeleton points (NeighborSkeletonIDs);
%               the IDs of neighboring tls points (NeighborTLSIDs); the sigma value of a skeleton point (Sigma); 
%               the average term used to calculate the coordinates of a skeleton (AverageTerm); and the repulsion term used to calculate the coordinates of a skeleton (RepulsionTerm)
% tlsPoints    The struct that records the attributes of TLS points.
%              The attributes include the ID of a TLS point (ID); the coordinates of a TLS point (P);
%              the IDs of neighboring tls points (NeighborTLSIDs); the coordinates of neighboring tls points (NeighborTLSPoints);
%              the density around a TLS point (Density)
% ------------------------------------------------------------------------------

function [seedPoints,tlsPoints] = L1TreeExtractSeedPoints(tlsPoints,paths,parameters)

if parameters.leafOn
    delIDs = find(tlsPoints(:,4) == 0);    % Delete TLS point with a frequency of 0
    if ~isempty(delIDs)
        tlsPoints(delIDs,:) = [];
    end
    tlsPoints = sortrows(tlsPoints,4,'descend');
    skeletonPointsNum = fix(size(tlsPoints,1)*parameters.skeletonProp);
    skPoints = tlsPoints(1:skeletonPointsNum,:);
    skNeighborIDs = rangesearch(skPoints(:,1:3),skPoints(:,1:3),parameters.nearThreshold*2,'Distance','euclidean','NSMethod','kdtree');
    delIDs = [];
    for i = 1:1:size(skPoints,1)
        if ismember(i,delIDs)
            continue;
        end
        delIDs = [delIDs,setdiff(skNeighborIDs{i},i)];
    end
    delIDs = sort(unique(delIDs));
    skPoints(delIDs,:) = [];
    skeletonPointsNum = size(skPoints,1);
    isTLS = false(size(tlsPoints,1),1);
    for i = 1:1:skeletonPointsNum
        tempPoint = skPoints(i,1:3);
        tempDistanceThreshold = max(sqrt(skPoints(i,6)/skPoints(i,5))*parameters.baseRadius*2,parameters.initialSearchRange*2);
        tempDistances = sqrt(sum((tlsPoints(:,1:3) - tempPoint).^2,2));
        isTLS(tempDistances < tempDistanceThreshold) = true;
    end
    tlsPoints = tlsPoints(isTLS,:);
end

if ~parameters.leafOn
    tlsNeighborMatrix = cell(size(tlsPoints,1),1);
    for i = 1:1:size(tlsPoints,1)
        tempPoint = tlsPoints(i,1:3);
        tempDistanceThreshold = max(sqrt(tlsPoints(i,6)/tlsPoints(i,5))*parameters.baseRadius*2,parameters.initialSearchRange*2);
        tempDistances = sqrt(sum((tlsPoints(:,1:3) - tempPoint).^2,2));
        selectIDs = find(tempDistances < tempDistanceThreshold);
        tlsNeighborMatrix(i) = {selectIDs};
    end
    skeletonPointIDs = [];
    isConsider = false(1,size(tlsPoints,1));
    isConsider(end) = true;
    while any(~isConsider)
        selectIDs = find(~isConsider);
        [~,idx] = min(sqrt(tlsPoints(selectIDs,6)./tlsPoints(selectIDs,5)));
        currentPointID = selectIDs(idx);
        addSkeletonPointIDs = paths{currentPointID};
        if isempty(addSkeletonPointIDs)
            isConsider(currentPointID) = true;
        else
            addSkeletonPointIDs = addSkeletonPointIDs(~ismember(addSkeletonPointIDs,skeletonPointIDs));
            skeletonPointIDs = [skeletonPointIDs,addSkeletonPointIDs];
            for i = 1:1:length(addSkeletonPointIDs)
                isConsider(tlsNeighborMatrix{addSkeletonPointIDs(i)}) = true;
            end
        end
    end
    skeletonPointIDs = sort(unique(skeletonPointIDs));
    skPoints = tlsPoints(skeletonPointIDs,:);
    skNeighborIDs = rangesearch(skPoints(:,1:3),skPoints(:,1:3),parameters.nearThreshold*2,'Distance','euclidean','NSMethod','kdtree');
    delIDs = [];
    for i = 1:1:size(skPoints,1)
        if ismember(i,delIDs)
            continue;
        end
        delIDs = [delIDs,setdiff(skNeighborIDs{i},i)];
    end
    delIDs = sort(unique(delIDs));
    skPoints(delIDs,:) = [];
    skeletonPointsNum = size(skPoints,1);
    delIDs = find(tlsPoints(:,4) == 0);
    if ~isempty(delIDs)
        tlsPoints(delIDs,:) = [];
    end
end

seedPoints = struct(...
    'ID',num2cell(linspace(1,skeletonPointsNum,skeletonPointsNum),1),...
    'P',num2cell(skPoints(:,1:3),2)',...
    'PathLength',num2cell(skPoints(:,5)',1),...
    'Location',num2cell(skPoints(:,6)',1),...
    'IsFix',num2cell(false(skeletonPointsNum,1)',1),...
    'IsIgnore',num2cell(false(skeletonPointsNum,1)',1),...
    'IsBranch',num2cell(false(skeletonPointsNum,1)',1),...
    'NeighborSkeletonIDs',cell(1,skeletonPointsNum),...
    'NeighborTLSIDs',cell(1,skeletonPointsNum),...
    'Sigma',num2cell(zeros(skeletonPointsNum,1)',1),...
    'AverageTerm',num2cell(zeros(skeletonPointsNum,3),2)',...
    'AverageWeightSum',num2cell(zeros(skeletonPointsNum,1)',1),...
    'RepulsionTerm',num2cell(zeros(skeletonPointsNum,3),2)',...
    'RepulsionWeightSum',num2cell(zeros(skeletonPointsNum,1)',1)...
);

tempPoints = struct(...
    'ID',num2cell(linspace(1,size(tlsPoints,1),size(tlsPoints,1)),1),...
    'P',num2cell(tlsPoints(:,1:3),2)',...
    'NeighborTLSIDs',cell(1,size(tlsPoints,1)),...
    'NeighborTLSPoints',cell(1,size(tlsPoints,1)),...
    'Density',num2cell(zeros(1,size(tlsPoints,1)),1)...
);
tlsPoints = tempPoints;

end