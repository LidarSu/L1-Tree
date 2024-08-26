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
% L1TreePruneBranches.m    The function for pruning terminal branches
%
% Version 1.0
% Latest update     04 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% treeModel    The input 3D tree model
% tlsPoints    The coordinates and attributes of TLS points
% parameters    The parameter set
% globalBar    The global waiting bar
% 
% OUTPUTS:
% treeModel    The updated 3D tree model 
% tlsPoints    The updated TLS points
% ------------------------------------------------------------------------------

function [treeModel,tlsPoints] = L1TreePruneBranches(treeModel,tlsPoints,parameters,globalBar)

waitbar(1/5,globalBar,'Prune terminal branches');

%% Step 1: Obtain the coordinates of skeleton points
skeletonPoints = [];
for i = 1:1:size(treeModel,2)
    skeletonPoints = [skeletonPoints;[treeModel(i).SkeletonIDs',treeModel(i).Curve]];
end
skeletonPoints = unique(skeletonPoints,'rows','stable');
skeletonPoints = sortrows(skeletonPoints,1);

%% Step 2: Identify the TLS points associated with each branch
selectIDs = find(tlsPoints(:,4) ~= 0);
associatedTLSPoints = tlsPoints(selectIDs,1:3);
tlsPointsToBranches = ones(size(associatedTLSPoints,1),2)*9999;
for i = 1:1:size(treeModel,2)
    for j = 1:1:length(treeModel(i).SkeletonIDs)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - associatedTLSPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - associatedTLSPoints).^2,2));
        distanceToPoint_1 = sqrt(sum((pedalPoints - point_1).^2,2));
        distanceToPoint_2 = sqrt(sum((pedalPoints - point_2).^2,2));
        selectIDs_1 = find((pedalDistances < tlsPointsToBranches(:,2)) & (distanceToPoint_1 < distanceToPoint_2));
        selectIDs_2 = find((pedalDistances < tlsPointsToBranches(:,2)) & (distanceToPoint_2 <= distanceToPoint_1));
        tlsPointsToBranches(selectIDs_1,1) = treeModel(i).SkeletonIDs(j);
        tlsPointsToBranches(selectIDs_1,2) = pedalDistances(selectIDs_1);
        tlsPointsToBranches(selectIDs_2,1) = treeModel(i).SkeletonIDs(j+1);
        tlsPointsToBranches(selectIDs_2,2) = pedalDistances(selectIDs_2);
    end
end

%% Step 3: Prune each terminal branch
prunedBranches = [];
count = 1;
for i = 1:1:size(treeModel,2)
    if ~treeModel(i).IsTerminal || (length(treeModel(i).SkeletonIDs) <= parameters.branchSizeThreshold)
        continue;
    end
    segRadii = zeros(1,length(treeModel(i).SkeletonIDs));
    for j = 1:1:length(treeModel(i).SkeletonIDs)
        if j == length(treeModel(i).SkeletonIDs)
            tempDirection = treeModel(i).Curve(j,:) - treeModel(i).Curve(j-1,:);
        else
            tempDirection = treeModel(i).Curve(j+1,:) - treeModel(i).Curve(j,:);
        end
        tempDirection = tempDirection/sqrt(sum(tempDirection.^2));
        tempSkeletonPointID = treeModel(i).SkeletonIDs(j);
        tempSkeletonPoint = treeModel(i).Curve(j,:);
        selectIDs = find(tlsPointsToBranches(:,1) == tempSkeletonPointID);
        tempTLSPoints = associatedTLSPoints(selectIDs,:);
        distances_1 = sqrt(sum((tempTLSPoints - tempSkeletonPoint).^2,2));
        distances_2 = abs((tempTLSPoints - tempSkeletonPoint)*tempDirection');
        segRadii(j) = mean(sqrt(distances_1.^2 - distances_2.^2));
    end
    breakID = L1TreeMovingTtest(real(segRadii));
    if isempty(breakID) || (breakID < parameters.branchSizeThreshold)
        continue;
    end
    prunedBranches(count).Curve = treeModel(i).Curve(breakID-1:end,:);
    count = count + 1;
    treeModel(i).SkeletonIDs = treeModel(i).SkeletonIDs(1:breakID-1);
    treeModel(i).Curve = treeModel(i).Curve(1:breakID-1,:);
end

%% Step 4: Update the 3D tree model
rootID = treeModel(1).SkeletonIDs(1);
neighborMatrix = [];
for i = 1:1:size(treeModel,2)
    for j = 1:1:length(treeModel(i).SkeletonIDs)-1
        pointID_1 = treeModel(i).SkeletonIDs(j);
        pointID_1 = find(skeletonPoints(:,1) == pointID_1);
        pointID_2 = treeModel(i).SkeletonIDs(j+1);
        pointID_2 = find(skeletonPoints(:,1) == pointID_2);
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        distance = sqrt(sum((point_1 - point_2).^2));
        neighborMatrix = [neighborMatrix;[pointID_1,pointID_2,distance]];
    end
end
treeModel = L1TreeUpdate(neighborMatrix,skeletonPoints(:,2:4),rootID);

%% Step 5: Update the attributes of each TLS point
tlsPointsToBranches = ones(size(tlsPoints,1),2)*9999;
for i = 1:1:size(treeModel,2)
    for j = 1:1:length(treeModel(i).SkeletonIDs)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - tlsPoints(:,1:3);
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - tlsPoints(:,1:3)).^2,2));
        selectIDs = find(pedalDistances < tlsPointsToBranches(:,2));
        tlsPointsToBranches(selectIDs,1) = i;
        tlsPointsToBranches(selectIDs,2) = pedalDistances(selectIDs);
    end
end
if isempty(prunedBranches)
    tlsPointsToPrunedBranches = ones(size(tlsPoints,1),1)*9999;
else
    tlsPointsToPrunedBranches = ones(size(tlsPoints,1),1)*9999;
    for i = 1:1:size(prunedBranches,2)
        for j = 1:1:size(prunedBranches(i).Curve,1)-1
            point_1 = prunedBranches(i).Curve(j,:);
            point_2 = prunedBranches(i).Curve(j+1,:);
            P1P2 = point_1 - point_2;
            P2P1 = point_2 - point_1;
            P1P0 = point_1 - tlsPoints(:,1:3);
            tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
            tMatrix(tMatrix < 0) = 0;
            tMatrix(tMatrix > 1) = 1;
            pedalPoints = point_1 + tMatrix*P2P1;
            pedalDistances = sqrt(sum((pedalPoints - tlsPoints(:,1:3)).^2,2));
            selectIDs = find(pedalDistances < tlsPointsToPrunedBranches);
            tlsPointsToPrunedBranches(selectIDs) = pedalDistances(selectIDs);
        end
    end
end

selectIDs_1 = find(tlsPointsToPrunedBranches < tlsPointsToBranches(:,2));
selectIDs_2 = find(tlsPoints(:,4) == 0);
tlsPoints(:,4) = tlsPointsToBranches(:,1);
tlsPoints(selectIDs_1,4) = 0;
tlsPoints(selectIDs_2,4) = 0;

end