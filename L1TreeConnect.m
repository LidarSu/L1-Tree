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
% L1TreeConnect.m    The function for connecting the broken skeleton
%
% Version 1.0
% Latest update     31 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% skeletonBranches    The initial skeleton with multiple branch segments
% skeletonPoints    The input skeleton points
% tlsPoints    The input tlsPoints
% parameters    The parameter set
% globalBar    The global waiting bar
%
% OUTPUTS:
% treeModel    The struct of the 3D tree model, which records the IDs of the skeleton points that make up each branch (SkeletonIDs)
%              the coordinates of the skeleton points that make up each branch (Curve); the centrifugal code of each branch (CenCode)
%              the centrifugal order of each branch (CenOrder); the radius of each branch under the centrifugal coding strategy (CenRadius)
%              the length of each branch under the centrifugal coding strategy (CenLength); 
%              whether the radius value is credible under the centrifugal coding strategy (the term "incredible" means that the radius of a branch is greater than that of its parent; IsCenCredible);
%              whether a branch is a terminal one under the centrifufal coding strategy (IsTerminal);
%              Gaaliche’s code of each branch (GaCode); Gaaliche’s oder of each branch (GaOrder);
%              the radius of each branch under Gaaliche’s coding strategy (GaRadius); the length of each branch under Gaaliche’s coding strategy (GaLength);
%              and wheter the radius value is credible under Gaaliche’s coding strategy (IsGaCredible)
% tlsPoints    The coordinates and attributes of TLS points.
%              The first three columns record the x, y, and z coordinates of TLS points; and the last column records the ID of the branch to which each TLS point belongs
% ------------------------------------------------------------------------------

function [treeModel,tlsPoints] = L1TreeConnect(skeletonBranches,skeletonPoints,tlsPoints,parameters,globalBar)

waitbar(0,globalBar,'Connect the broken skeleton');

%% Step 1: Identify TLS points, fixed skeleton points, unfixed skeleton points, and end points of each branch segments
% Step 1.1: Identify TLS points
tempPoints = [tlsPoints(:).P];
tlsPoints = reshape(tempPoints,3,size(tlsPoints,2))';
% Step 1.2: Remove unreasonable branch segments
% Step 1.2.1: Calculate the distance from each TLS point to its nearest branch segment
tlsPointsToBranches = ones(size(tlsPoints,1),4)*9999;
for i = 1:1:size(skeletonBranches,2)
    for j = 1:1:size(skeletonBranches(i).Curve,1)-1
        point_1 = skeletonBranches(i).Curve(j,:);
        point_2 = skeletonBranches(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - tlsPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix > 1) = 1;
        tMatrix(tMatrix < 0) = 0;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
        selectIDs = find(pedalDistances < tlsPointsToBranches(:,4));
        tlsPointsToBranches(selectIDs,1) = i;
        tlsPointsToBranches(selectIDs,2) = j;
        if j == 1
            selectIDs_1 = find((tMatrix == 0) & (pedalDistances < tlsPointsToBranches(:,4)));
            selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < tlsPointsToBranches(:,4)));
            tlsPointsToBranches(selectIDs_1,3) = 1;
            tlsPointsToBranches(selectIDs_2,3) = 0;
        elseif j == size(skeletonBranches(i).Curve,1)-1
            selectIDs_1 = find((tMatrix == 1) & (pedalDistances < tlsPointsToBranches(:,4)));
            selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < tlsPointsToBranches(:,4)));
            tlsPointsToBranches(selectIDs_1,3) = 1;
            tlsPointsToBranches(selectIDs_2,3) = 0;
        else
            tlsPointsToBranches(selectIDs,3) = 0;
        end
        tlsPointsToBranches(selectIDs,4) = pedalDistances(selectIDs);
    end
end
% Step 1.2.2: Remove the first type of branches. These branches do not pass through the dense point cloud
delIDs_1 = [];
for i = 1:1:size(skeletonBranches,2)
    totalLength = 0;
    maxSegmentLength = 0;
    maxSegmentID = 0;
    for j = 1:1:size(skeletonBranches(i).Curve,1)-1
        point_1 = skeletonBranches(i).Curve(j,:);
        point_2 = skeletonBranches(i).Curve(j+1,:);
        segmentLength = sqrt(sum((point_1 - point_2).^2));
        if segmentLength > maxSegmentLength
            maxSegmentLength = segmentLength;
            maxSegmentLength = j;
        end
        totalLength = totalLength + segmentLength;
    end
    if (totalLength <= parameters.minDistanceThreshold) || (maxSegmentLength <= parameters.windowSize)    % Here used to be "totalLength <= parameters.minTestLinkLength = 0.15"
        continue;
    end
    branchRelatedTLSPointIDs = find((tlsPointsToBranches(:,1) == i) & (tlsPointsToBranches(:,3) == 0));
    tempSearchRange = quantile(tlsPointsToBranches(branchRelatedTLSPointIDs,4),0.9);
    branchTLSPointNum = length(branchRelatedTLSPointIDs);
    segmentTLSPointNum = find((tlsPointsToBranches(:,1) == i) & (tlsPointsToBranches(:,2) == maxSegmentID) & (tlsPointsToBranches(:,3) == 0) & (tlsPointsToBranches(:,4) < tempSearchRange));
    if segmentTLSPointNum < branchTLSPointNum/totalLength*maxSegmentLength*parameters.densityFactor
        delIDs_1 = [delIDs_1,i];
    end
end
% Step 1.2.3: Remove the second type of branches. These branches are not in the center of the point cloud
delIDs_2 = [];
for i = 1:1:size(skeletonBranches,2)
    branchRelatedTLSPointIDs = find((tlsPointsToBranches(:,1) == i) & (tlsPointsToBranches(:,3) == 0) & (tlsPointsToBranches(:,4) > 0.000001));
    if isempty(branchRelatedTLSPointIDs)
        delIDs_2 = [delIDs_2,i];
        continue;
    end
    branchRelatedTLSPoints = tlsPoints(branchRelatedTLSPointIDs,:);
    pedalInfo = ones(size(branchRelatedTLSPoints,1),4)*9999;
    for j = 1:1:size(skeletonBranches(i).Curve,1)-1
        point_1 = skeletonBranches(i).Curve(j,:);
        point_2 = skeletonBranches(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - branchRelatedTLSPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix > 1) = 1;
        tMatrix(tMatrix < 0) = 0;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - branchRelatedTLSPoints).^2,2));
        selectIDs = find(pedalDistances < pedalInfo(:,4));
        pedalInfo(selectIDs,1:3) = pedalPoints(selectIDs,1:3);
        pedalInfo(selectIDs,4) = pedalDistances(selectIDs);
    end
    branchDirection = skeletonBranches(i).Curve(end,:) - skeletonBranches(i).Curve(1,:);
    branchDirection = branchDirection/sqrt(sum(branchDirection.^2));
    oldAxisX = [1,0,0];
    oldAxisY = [0,1,0];
    oldAxisZ = [0,0,1];
    A = [oldAxisX;oldAxisY;oldAxisZ];
    oldCoordinates = skeletonBranches(i).Curve(1,:) + branchRelatedTLSPoints - pedalInfo(:,1:3);
    newAxisX = branchRelatedTLSPoints(1,:) - pedalInfo(1,1:3);
    newAxisX = newAxisX/sqrt(sum(newAxisX.^2));
    newAxisZ = branchDirection;
    newAxisY = cross(newAxisX,newAxisZ);
    newAxisY = newAxisY/sqrt(sum(newAxisY.^2));
    B = [newAxisX;newAxisY;newAxisZ];
    newCoordinates = (oldCoordinates*A - skeletonBranches(i).Curve(1,:))*inv(B);
    newCoordinates = newCoordinates(:,1:2);
    distances = sqrt(sum(newCoordinates.^2,2));
    newCoordinates(:,1) = newCoordinates(:,1)./distances;
    newCoordinates(:,2) = newCoordinates(:,2)./distances;
    angles = ones(1,size(newCoordinates,1))*9999;
    for j = 1:1:size(newCoordinates,1)
        angle = acos(newCoordinates(j,1))*180/pi;
        if newCoordinates(j,2) < 0
            angle = 360 - angle;
        end
        angles(j) = angle;
    end
    angles = sort(angles);
    angleDiff = ones(1,size(newCoordinates,1))*9999;
    for j = 1:1:size(newCoordinates,1)
        if j == 1
            angleDiff(j) = angles(1) + 360 - angles(end);
        else
            angleDiff(j) = angles(j) - angles(j-1);
        end
    end
    angleDiff = sort(angleDiff,'descend');
    if(size(newCoordinates,1) > 1) && (angleDiff(2)*sqrt(size(newCoordinates,1)) > 720)
        delIDs_2 = [delIDs_2,i];
    end
end
% Step 1.2.4: Remove the third type of branches. These branches are within the inflence of other branches
branchLengthes = ones(1,size(skeletonBranches,2))*9999;
for i = 1:1:size(skeletonBranches,2)
    totalLength = 0;
    for j = 1:1:size(skeletonBranches(i).Curve,1)-1
        point_1 = skeletonBranches(i).Curve(j,:);
        point_2 = skeletonBranches(i).Curve(j+1,:);
        segmentLength = sqrt(sum((point_1 - point_2).^2));
        totalLength = totalLength + segmentLength;
    end
    branchLengthes(i) = totalLength;
end
[~,idx] = sort(branchLengthes,'descend');
delIDs_3 = [];
for i = 1:1:length(idx)
    if ismember(idx(i),delIDs_1) || ismember(idx(i),delIDs_2) || ismember(idx(i),delIDs_3)
        continue;
    end
    currentBranch = skeletonBranches(idx(i));
    branchRelatedTLSPointIDs = find((tlsPointsToBranches(:,1) == idx(i)) & (tlsPointsToBranches(:,2) == 0));
    tempSearchRange = quantile(tlsPointsToBranches(branchRelatedTLSPointIDs,4),0.9);
    for j = 1:1:size(skeletonBranches,2)
        if ismember(idx(j),delIDs_1) || ismember(idx(j),delIDs_2) || ismember(idx(j),delIDs_3) || (j == idx(i))
            continue;
        end
        consideredBranchPoints = skeletonBranches(j).Curve;
        consideredToCurrentBranch = ones(size(consideredBranchPoints,1),2)*9999;
        for k = 1:1:size(currentBranch.Curve,1)-1
            point_1 = currentBranch.Curve(k,:);
            point_2 = currentBranch.Curve(k+1,:);
            P1P2 = point_1 - point_2;
            P2P1 = point_2 - point_1;
            P1P0 = point_1 - consideredBranchPoints;
            tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
            tMatrix(tMatrix > 1) = 1;
            tMatrix(tMatrix < 0) = 0;
            pedalPoints = point_1 + tMatrix*P2P1;
            pedalDistances = sqrt(sum((pedalPoints - consideredBranchPoints).^2,2));
            selectIDs = find(pedalDistances < consideredToCurrentBranch(:,2));
            if k == 1
                selectIDs_1 = find((tMatrix == 0) & (pedalDistances < consideredToCurrentBranch(:,2)));
                selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < consideredToCurrentBranch(:,2)));
                consideredToCurrentBranch(selectIDs_1,1) = 1;
                consideredToCurrentBranch(selectIDs_2,1) = 0;
            elseif k == size(currentBranch.Curve,1)-1
                selectIDs_1 = find((tMatrix == 1) & (pedalDistances < consideredToCurrentBranch(:,2)));
                selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < consideredToCurrentBranch(:,2)));
                consideredToCurrentBranch(selectIDs_1,1) = 1;
                consideredToCurrentBranch(selectIDs_2,1) = 0;
            else
                consideredToCurrentBranch(selectIDs,1) = 0;
            end
            consideredToCurrentBranch(selectIDs,2) = pedalDistances(selectIDs);
        end
        selectIDs = find(consideredToCurrentBranch(:,1) == 0);
        if isempty(selectIDs)
            continue;
        end
        if (mean(consideredToCurrentBranch(selectIDs,2)) < tempSearchRange) && (mean(consideredToCurrentBranch(selectIDs,2)) < parameters.baseRadius) && (length(selectIDs) > size(consideredBranchPoints,1)/2)
            delIDs_3 = [delIDs_3,j];
        end
    end
end
% Step 1.3: Identify fixed and unfixed skeleton points. 
% The skeleton points on branches with delIDs_1 are considered as unfixed and will be involved in subsequent calculations; the skeleton points on branches with delIDs_2 and delIDs_3 are considered fixed and will not be involved in subsequent calculations
% Step 1.3.1: Identify unfixed skeleton points
selectIDs = find(~[skeletonPoints(:).IsIgnore] & ~[skeletonPoints(:).IsFix]);
unfixedPoints = [skeletonPoints(selectIDs).P];
unfixedPoints = reshape(unfixedPoints,3,length(selectIDs))';
for i = 1:1:length(delIDs_1)
    unfixedPoints = [unfixedPoints;skeletonBranches(delIDs_1(i)).Curve];
end
unfixedPoints = unique(unfixedPoints,'rows','stable');
delIDs = unique([delIDs_1,delIDs_2,delIDs_3]);
skeletonBranches(delIDs) = [];
% Step 1.3.2: Identify fixed skeleton points
occurrenceOrder = [];
for i = 1:1:size(skeletonBranches,2)
    for j = 1:1:length(skeletonBranches(i).SkeletonIDs)
        if ~ismember(skeletonBranches(i).SkeletonIDs(j),occurrenceOrder)
            occurrenceOrder = [occurrenceOrder,skeletonBranches(i).SkeletonIDs(j)];
        end
    end
end
for i = 1:1:size(skeletonBranches,2)
    newSkeletonIDs = [];
    for j = 1:1:length(skeletonBranches(i).SkeletonIDs)
        newSkeletonIDs = [newSkeletonIDs,find(occurrenceOrder == skeletonBranches(i).SkeletonIDs(j))];
    end
    skeletonBranches(i).SkeletonIDs = newSkeletonIDs;
end
tempPoints = [skeletonPoints(occurrenceOrder).P];
skeletonPoints = reshape(tempPoints,3,length(occurrenceOrder))';
allSkeletonPoints = [skeletonPoints;unfixedPoints];
% Step 1.4: Identify end points of branch segments
endPoints = [];
for i = 1:1:size(skeletonBranches,2)
    if length(skeletonBranches(i).SkeletonIDs) < 4
        endPoint.BranchID = i;
        endPoint.ID = skeletonBranches(i).SkeletonIDs(1);
        endPoint.P = allSkeletonPoints(endPoint.ID,:);
        endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(1),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end),:);
        endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
        endPoints = [endPoints,endPoint];
        endPoint.BranchID = i;
        endPoint.ID = skeletonBranches(i).SkeletonIDs(end);
        endPoint.P = allSkeletonPoints(endPoint.ID,:);
        endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(1),:);
        endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
        endPoints = [endPoints,endPoint];
    else
        point_1 = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(1),:);
        point_2 = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(2),:);
        point_3 = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end-1),:);
        point_4 = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end),:);
        direction_1 = (point_2 - point_1)/sqrt(sum((point_2 - point_1).^2));
        direction_2 = (point_4 - point_3)/sqrt(sum((point_4 - point_3).^2));
        angle = acos(direction_1*direction_2')*180/pi;
        if angle < 90
            endPoint.BranchID = i;
            endPoint.ID = skeletonBranches(i).SkeletonIDs(1);
            endPoint.P = allSkeletonPoints(endPoint.ID,:);
            endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(1),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end),:);
            endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
            endPoints = [endPoints,endPoint];
            endPoint.BranchID = i;
            endPoint.ID = skeletonBranches(i).SkeletonIDs(end);
            endPoint.P = allSkeletonPoints(endPoint.ID,:);
            endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(1),:);
            endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
            endPoints = [endPoints,endPoint];
        else
            angles = [];
            for j = 2:1:length(skeletonBranches(i).SkeletonIDs)-2
                point_2 = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(j),:);
                point_3 = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(j+1),:);
                direction_1 = (point_2 - point_1)/sqrt(sum((point_2 - point_1).^2));
                direction_2 = (point_4 - point_3)/sqrt(sum((point_4 - point_3).^2));
                angle = acos(direction_1*direction_2')*180/pi;
                angles = [angles,angle];
            end
            [~,maxID] = max(angles);
            endPoint.BranchID = i;
            endPoint.ID = skeletonBranches(i).SkeletonIDs(1);
            endPoint.P = allSkeletonPoints(endPoint.ID,:);
            endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(1),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(maxID+1),:);
            endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
            endPoints = [endPoints,endPoint];
            endPoint.BranchID = i;
            endPoint.ID = skeletonBranches(i).SkeletonIDs(maxID+1);
            endPoint.P = allSkeletonPoints(endPoint.ID,:);
            endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(maxID+1),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(1),:);
            endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
            endPoints = [endPoints,endPoint];
            endPoint.BranchID = i;
            endPoint.ID = skeletonBranches(i).SkeletonIDs(maxID+2);
            endPoint.P = allSkeletonPoints(endPoint.ID,:);
            endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(maxID+2),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end),:);
            endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
            endPoints = [endPoints,endPoint];
            endPoint.BranchID = i;
            endPoint.ID = skeletonBranches(i).SkeletonIDs(end);
            endPoint.P = allSkeletonPoints(endPoint.ID,:);
            endPoint.Direction = allSkeletonPoints(skeletonBranches(i).SkeletonIDs(end),:) - allSkeletonPoints(skeletonBranches(i).SkeletonIDs(maxID+2),:);
            endPoint.Direction = endPoint.Direction/sqrt(sum(endPoint.Direction.^2));
            endPoints = [endPoints,endPoint];
        end
    end
end
% Step 1.5: Recalculate the distances from TLS points to branches
tlsPointsToBranches = ones(size(tlsPoints,1),2)*9999;
for i = 1:1:size(skeletonBranches,2)
    for j = 1:1:size(skeletonBranches(i).Curve,1)-1
        point_1 = skeletonBranches(i).Curve(j,:);
        point_2 = skeletonBranches(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - tlsPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
        selectIDs = find(pedalDistances < tlsPointsToBranches(:,2));
        if j  == 1
            selectIDs_1 = find((tMatrix == 0) & (pedalDistances < tlsPointsToBranches(:,2)));
            selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < tlsPointsToBranches(:,2)));
            tlsPointsToBranches(selectIDs_1,1) = 1;
            tlsPointsToBranches(selectIDs_2,1) = 0;
        elseif j == size(skeletonBranches(i).Curve,1)-1
            selectIDs_1 = find((tMatrix == 1) & (pedalDistances < tlsPointsToBranches(:,2)));
            selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < tlsPointsToBranches(:,2)));
            tlsPointsToBranches(selectIDs_1,1) = 1;
            tlsPointsToBranches(selectIDs_2,1) = 0;
        else
            tlsPointsToBranches(selectIDs,1) = 0;
        end
        tlsPointsToBranches(selectIDs,2) = pedalDistances(selectIDs);
    end
end

%% Step 2: Generate the neighboring matrix
% Step 2.1: Generate the neighboring matrix between fixed points
neighborMatrix = [];
for i = 1:1:size(skeletonBranches,2)
    for j = 1:1:length(skeletonBranches(i).SkeletonIDs)-1
        pointID_1 = skeletonBranches(i).SkeletonIDs(j);
        pointID_2 = skeletonBranches(i).SkeletonIDs(j+1);
        point_1 = allSkeletonPoints(pointID_1,:);
        point_2 = allSkeletonPoints(pointID_2,:);
        distance = sqrt(sum((point_1 - point_2).^2));
        neighborMatrix = [neighborMatrix;[pointID_1,pointID_2,distance]];
    end
end
ignoreDirections = [];
for i = 1:1:size(endPoints,2)-1
    for j = (i+1):1:size(endPoints,2)
        pointID_1 = endPoints(i).ID;
        pointID_2 = endPoints(j).ID;
        point_1 = allSkeletonPoints(pointID_1,:);
        point_2 = allSkeletonPoints(pointID_2,:);
        distance = sqrt(sum((point_1 - point_2).^2));
        if distance <= parameters.minDistanceThreshold
            neighborMatrix = [neighborMatrix;[pointID_1,pointID_2,distance]];
            ignoreDirections = [ignoreDirections;[pointID_1,pointID_2];[pointID_2,pointID_1]];
        end
    end
end
% Step 2.2: Generate the neighboring matrix related to unfixed skeleton points
unfixedNeighborMatrix = [];
for i = 1:1:size(unfixedPoints,1)-1
    for j = (i+1):1:size(unfixedPoints,1)
        point_1 = unfixedPoints(i,:);
        point_2 = unfixedPoints(j,:);
        distance = sqrt(sum((point_1 - point_2).^2));
        if distance <= parameters.minDistanceThreshold
            pointID_1 = find((allSkeletonPoints(:,1) == point_1(1)) & (allSkeletonPoints(:,2) == point_1(2)) & (allSkeletonPoints(:,3) == point_1(3)));
            pointID_2 = find((allSkeletonPoints(:,1) == point_2(1)) & (allSkeletonPoints(:,2) == point_2(2)) & (allSkeletonPoints(:,3) == point_2(3)));
            unfixedNeighborMatrix = [unfixedNeighborMatrix;[pointID_1,pointID_2,distance]];
        end
    end
end
for i = 1:1:size(unfixedPoints,1)
    for j = 1:1:size(skeletonPoints,1)
        point_1 = unfixedPoints(i,:);
        point_2 = skeletonPoints(j,:);
        distance = sqrt(sum((point_1 - point_2).^2));
        if distance <= parameters.minDistanceThreshold
            pointID_1 = find((allSkeletonPoints(:,1) == point_1(1)) & (allSkeletonPoints(:,2) == point_1(2)) & (allSkeletonPoints(:,3) == point_1(3)));
            pointID_2 = find((allSkeletonPoints(:,1) == point_2(1)) & (allSkeletonPoints(:,2) == point_2(2)) & (allSkeletonPoints(:,3) == point_2(3)));
            unfixedNeighborMatrix = [unfixedNeighborMatrix;[pointID_1,pointID_2,distance]];
        end
    end
end
potentialConnections = graph(unfixedNeighborMatrix(:,1),unfixedNeighborMatrix(:,2),unfixedNeighborMatrix(:,3));

%% Step 3: Generate connected components
treeGraph = graph(neighborMatrix(:,1),neighborMatrix(:,2),neighborMatrix(:,3));
[~,pred] = minspantree(treeGraph,'Method','dense','Type','forest');
rootIDs = find(pred == 0);
count = 1;
for i = 1:1:length(rootIDs)
    % Step 3.1: Identify the root point of each connected component
    rootID = rootIDs(i);
    nextIDs = find(pred == rootID);
    if isempty(nextIDs)
        continue;
    end
    tempComponent = [];
    if length(nextIDs) == 1
        nextID = nextIDs(1);
        currentCenCode = 'T';
        tempComponent = L1TreeGetCenCode(tempComponent,pred,currentCenCode,rootID,nextID);
    else
        for j = 1:1:length(nextIDs)
            nextID = nextIDs(j);
            currentCenCode = ['T_B',num2str(j)];
            tempComponentAdd = [];
            tempComponentAdd = L1TreeGetCenCode(tempComponentAdd,pred,currentCenCode,rootID,nextID);
            tempComponent = [tempComponent,tempComponentAdd];
        end
    end
    potentialRootIDs = [];
    for j = 1:1:size(tempComponent,2)
        if tempComponent(j).CenOrder == min([tempComponent(:).CenOrder])
            if length(find([tempComponent(:).SkeletonIDs] == tempComponent(j).SkeletonIDs(1))) == 1
                potentialRootIDs = [potentialRootIDs,tempComponent(j).SkeletonIDs(1)];
            end
            if length(find([tempComponent(:).SkeletonIDs] == tempComponent(j).SkeletonIDs(end))) == 1
                potentialRootIDs = [potentialRootIDs,tempComponent(j).SkeletonIDs(end)];
            end
        end
        isEnd = true;
        for k = 1:1:size(tempComponent,2)
            if contains(tempComponent(k).CenCode,tempComponent(j).CenCode) && (tempComponent(k).CenOrder == tempComponent(j).CenOrder+1)
                isEnd = false;
            end
        end
        if isEnd
            if length(find([tempComponent(:).SkeletonIDs] == tempComponent(j).SkeletonIDs(1))) == 1
                potentialRootIDs = [potentialRootIDs,tempComponent(j).SkeletonIDs(1)];
            end
            if length(find([tempComponent(:).SkeletonIDs] == tempComponent(j).SkeletonIDs(end))) == 1
                potentialRootIDs = [potentialRootIDs,tempComponent(j).SkeletonIDs(end)];
            end
        end
    end
    potentialRootIDs = unique(potentialRootIDs);
    potentialRootPoints = allSkeletonPoints(potentialRootIDs,:);
    selectID = find(potentialRootPoints(:,3) == min(potentialRootPoints(:,3)),1);
    rootID = potentialRootIDs(selectID);
    % Step 3.2: Recode each branch of the connected component
    [~,newPred] = minspantree(treeGraph,'Method','dense','Type','tree','Root',rootID);
    nextIDs = find(newPred == rootID);
    tempComponent = [];
    if length(nextIDs) == 1
        nextID = nextIDs(1);
        currentCenCode = 'T';
        tempComponent = L1TreeGetCenCode(tempComponent,newPred,currentCenCode,rootID,nextID);
    else
        for j = 1:1:length(nextIDs)
            nextID = nextIDs(j);
            currentCenCode = ['T_B',num2str(j)];
            tempComponentAdd = [];
            tempComponentAdd = L1TreeGetCenCode(tempComponentAdd,newPred,currentCenCode,rootID,nextID);
            tempComponent = [tempComponent,tempComponentAdd];
        end
    end
    Components(count).ID = i;
    Components(count).Component = tempComponent;
    count = count + 1;
end
% Step 3.3: Sort connected components
ComponentLocations = [];
for i = 1:1:size(Components,2)
    tempComponent = Components(i).Component;
    tempLocation = 9999;
    for j = 1:1:size(tempComponent,2)
        if min(allSkeletonPoints(tempComponent(j).SkeletonIDs,3)) < tempLocation
            tempLocation = min(allSkeletonPoints(tempComponent(j).SkeletonIDs,3));
        end
    end
    ComponentLocations = [ComponentLocations,tempLocation];
end
[~,idx] = sort(ComponentLocations);
Components = Components(idx);
for i = 1:1:size(Components,2)
    Components(i).ID = i;
end
% Step 3.4: Calculate the distances from TLS points to connected components
tlsPointsToComponents = ones(size(tlsPoints,1),4)*9999;
for i = 1:1:size(Components,2)
    tempComponentID = Components(i).ID;
    tempComponent = Components(i).Component;
    for j = 1:1:size(tempComponent,2)
        for k = 1:1:length(tempComponent(j).SkeletonIDs)-1
            point_1 = allSkeletonPoints(tempComponent(j).SkeletonIDs(k),:);
            point_2 = allSkeletonPoints(tempComponent(j).SkeletonIDs(k+1),:);
            P1P2 = point_1 - point_2;
            P2P1 = point_2 - point_1;
            P1P0 = point_1 - tlsPoints;
            tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
            tMatrix(tMatrix < 0) = 0;
            tMatrix(tMatrix > 1) = 1;
            pedalPoints = point_1 + tMatrix*P2P1;
            pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
            selectIDs = find(pedalDistances < tlsPointsToComponents(:,4));
            tlsPointsToComponents(selectIDs,1) = tempComponentID;
            tlsPointsToComponents(selectIDs,2) = j;
            tlsPointsToComponents(selectIDs,4) = pedalDistances(selectIDs);
        end
    end
end
tlsPointsToComponents(:,3) = tlsPointsToBranches(:,1);

%% Step 4: Connect components that are far apart
oriComponentNum = size(Components,2);
targetComponent = Components(1).Component;
lastComponent = Components(1).Component;
lastComponentID = Components(1).ID;
Components(1) = [];
connectionStrategies = [];
allInsertPointIDs = [];
while true
    while ~isempty(Components)
        % Step 4.1: Get the information about the connected component
        lastComponentPointIDs = [];
        lastComponentPoints = [];
        lastRelatedTLSPointIDs = [];
        lastSearchRanges = [];
        for i = 1:1:size(lastComponent,2)
            if i == 1
                tempIDs = lastComponent(i).SkeletonIDs;
            else
                tempIDs = lastComponent(i).SkeletonIDs(2:end);
            end
            tempPoints = allSkeletonPoints(tempIDs,:);
            lastComponentPointIDs = [lastComponentPointIDs,tempIDs];
            lastComponentPoints = [lastComponentPoints;tempPoints];
            selectIDs = find((tlsPointsToComponents(:,1) == lastComponentID) & (tlsPointsToComponents(:,2) == i) & (tlsPointsToComponents(:,3) == 0));
            if ~isempty(selectIDs)
                searchRange = quantile(tlsPointsToComponents(selectIDs,4),0.9);
            else
                selectIDs = find((tlsPointsToComponents(:,1) == lastComponentID) & (tlsPointsToComponents(:,3) == 0));
                searchRange = quantile(tlsPointsToComponents(selectIDs,4),0.9);
            end
            lastRelatedTLSPointIDs = [lastRelatedTLSPointIDs;rangesearch(tlsPoints,tempPoints,searchRange)];
            lastSearchRanges = [lastSearchRanges,repelem(searchRange,size(tempPoints,1))];
        end
        % Step 4.2: Iterative over all possible connections, which are from the end points of the current component and all skeleton points of the last component
        for i = 1:1:size(Components,2)
            currentComponent = Components(i).Component;
            currentComponentID = Components(i).ID;
            currentComponentEndPoints = [];
            for j = 1:1:size(currentComponent,2)
                for k = 1:1:length(currentComponent(j).SkeletonIDs)
                    selectID = find([endPoints(:).ID] == currentComponent(j).SkeletonIDs(k));
                    if ~isempty(selectID)
                        tempPointID = endPoints(selectID).ID;
                        tempPoint = endPoints(selectID).P;
                        tempDirection = endPoints(selectID).Direction;
                        currentComponentEndPoints = [currentComponentEndPoints;[tempPointID,tempPoint,tempDirection]];
                    end
                end
            end
            for j = 1:1:size(currentComponentEndPoints,1)
                currentComponentEndPointID = currentComponentEndPoints(j,1);
                currentComponentEndPoint = currentComponentEndPoints(j,2:4);
                startDirection = currentComponentEndPoints(j,5:7);
                for k = 1:1:size(lastComponentPoints,1)
                    lastComponentPointID = lastComponentPointIDs(k);
                    lastComponentPoint = lastComponentPoints(k,:);
                    connectionParameters.CurrentComponentID = currentComponentID;
                    connectionParameters.CurrentComponent = currentComponent;
                    connectionParameters.StartDirection = startDirection;
                    connectionParameters.LastAddTLSPointIDs = lastRelatedTLSPointIDs{k}';
                    connectionParameters.LastSearchRange = lastSearchRanges(k);
                    connectionLength = sqrt(sum((lastComponentPoint - currentComponentEndPoint).^2));
                    if connectionLength > parameters.maxDistanceThreshold
                        continue;
                    end
                    % Determine whether a direct connection is credible
                    potentialConnection_1 = [currentComponentEndPointID,lastComponentPointID];
                    [isCredible,connectionStrategy] = L1TreeIsConnectionCredible(potentialConnection_1,targetComponent,allSkeletonPoints,tlsPoints,ignoreDirections,tlsPointsToComponents,connectionParameters,parameters);
                    if isCredible
                        connectionStrategies = [connectionStrategies,connectionStrategy];
                    end
                    % Determine whether the connection passing through unfixed skeleton points is credible
                    potentialGraphSize = max([unfixedNeighborMatrix(:,1)',unfixedNeighborMatrix(:,2)']);
                    if (currentComponentEndPointID > potentialGraphSize) || (lastComponentPointID > potentialGraphSize)
                        potentialConnection_2 = [];
                    else
                        potentialConnection_2 = shortestpath(potentialConnections,currentComponentEndPointID,lastComponentPointID);
                    end
                    if isempty(potentialConnection_2)
                        continue;
                    end
                    [isCredible,connectionStrategy] = L1TreeIsConnectionCredible(potentialConnection_2,targetComponent,allSkeletonPoints,tlsPoints,ignoreDirections,tlsPointsToComponents,connectionParameters,parameters);
                    if isCredible
                        connectionStrategies = [connectionStrategies,connectionStrategy];
                    end
                end
            end
        end
        if isempty(connectionStrategies)
            break;
        else
            parameters.maxDistanceThreshold = 0.8;
        end
        % Step 4.3: Identify the optimal connection strategy
        componentIDs = unique([connectionStrategies(:).CurrentComponentID]);
        optConnectionStrategies = [];
        for i = 1:1:length(componentIDs)
            selectIDs = find([connectionStrategies(:).CurrentComponentID] == componentIDs(i));
            tempConnectionStrategies = connectionStrategies(selectIDs);
            [~,maxID] = max([tempConnectionStrategies(:).RelatedNum]./[tempConnectionStrategies(:).Length]./[tempConnectionStrategies(:).Distance]);
            optConnectionStrategies = [optConnectionStrategies,tempConnectionStrategies(maxID)];
        end
        [~,maxID] = max([optConnectionStrategies(:).RelatedNum]./[optConnectionStrategies(:).Length]./[optConnectionStrategies(:).Distance]);
        optConnectionStrategy = optConnectionStrategies(maxID);
        % Step 4.4: Apply the optimal connection
        targetNeighborMatrix = [];
        lastNeighborMatrix = [];
        addSkeletonPoints = [];
        for i = 1:1:size(targetComponent,2)
            for j = 1:1:length(targetComponent(i).SkeletonIDs)-1
                pointID_1 = targetComponent(i).SkeletonIDs(j);
                pointID_2 = targetComponent(i).SkeletonIDs(j+1);
                point_1 = allSkeletonPoints(pointID_1,:);
                point_2 = allSkeletonPoints(pointID_2,:);
                distance = sqrt(sum((point_1 - point_2).^2));
                targetNeighborMatrix = [targetNeighborMatrix;[pointID_1,pointID_2,distance]];
            end
        end
        for i = 1:1:length(optConnectionStrategy.Connection)-1
            pointID_1 = optConnectionStrategy.Connection(i);
            pointID_2 = optConnectionStrategy.Connection(i+1);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            distance = sqrt(sum((point_1 - point_2).^2));
            if distance > parameters.resolution
                sampleNum = floor(distance/parameters.resolution);
                insertPoints = [];
                for j = 1:1:sampleNum
                    tempPoint = point_1 + (point_2 - point_1)*j/(sampleNum+1);
                    insertPoints = [insertPoints;tempPoint];
                end
                insertPoints = [point_1;insertPoints;point_2];
                for j = 1:1:size(insertPoints,1)-1
                    point_1 = insertPoints(j,:);
                    point_2 = insertPoints(j+1,:);
                    pointID_1 = find((allSkeletonPoints(:,1) == point_1(1)) & (allSkeletonPoints(:,2) == point_1(2)) & (allSkeletonPoints(:,3) == point_1(3)));
                    pointID_2 = find((allSkeletonPoints(:,1) == point_2(1)) & (allSkeletonPoints(:,2) == point_2(2)) & (allSkeletonPoints(:,3) == point_2(3)));
                    if isempty(pointID_1)
                        allSkeletonPoints = [allSkeletonPoints;point_1];
                        addSkeletonPoints = [addSkeletonPoints;point_1];
                        pointID_1 = size(allSkeletonPoints,1);
                        allInsertPointIDs = [allInsertPointIDs,pointID_1];
                    end
                    if isempty(pointID_2)
                        allSkeletonPoints = [allSkeletonPoints;point_2];
                        addSkeletonPoints = [addSkeletonPoints;point_2];
                        pointID_2 = size(allSkeletonPoints,1);
                        allInsertPointIDs = [allInsertPointIDs,pointID_2];
                    end
                    distance = sqrt(sum((point_1 - point_2).^2));
                    targetNeighborMatrix = [targetNeighborMatrix;[pointID_1,pointID_2,distance]];
                    lastNeighborMatrix = [lastNeighborMatrix;[pointID_1,pointID_2,distance]];
                end
            else
                targetNeighborMatrix = [targetNeighborMatrix;[pointID_1,pointID_2,distance]];
                lastNeighborMatrix = [lastNeighborMatrix;[pointID_1,pointID_2,distance]];
            end
        end
        selectID = find([Components(:).ID] == optConnectionStrategy.CurrentComponentID);
        optConnectionComponent = Components(selectID).Component;
        for i = 1:1:size(optConnectionComponent,2)
            for j = 1:1:length(optConnectionComponent(i).SkeletonIDs)-1
                pointID_1 = optConnectionComponent(i).SkeletonIDs(j);
                pointID_2 = optConnectionComponent(i).SkeletonIDs(j+1);
                point_1 = allSkeletonPoints(pointID_1,:);
                point_2 = allSkeletonPoints(pointID_2,:);
                distance = sqrt(sum((point_1 - point_2).^2));
                targetNeighborMatrix = [targetNeighborMatrix;[pointID_1,pointID_2,distance]];
                lastNeighborMatrix = [lastNeighborMatrix;[pointID_1,pointID_2,distance]];
            end
        end
        uniquePointIDs = unique([targetNeighborMatrix(:,1)',targetNeighborMatrix(:,2)']);
        uniquePoints = allSkeletonPoints(uniquePointIDs,:);
        [~,minID] = min(uniquePoints(:,3));
        rootID = uniquePointIDs(minID);
        treeGraph = graph(targetNeighborMatrix(:,1),targetNeighborMatrix(:,2),targetNeighborMatrix(:,3));
        [~,pred] = minspantree(treeGraph,'Method','dense','Type','tree','Root',rootID);
        nextID = find(pred == rootID);
        targetComponent = [];
        currentCenCode = 'T';
        targetComponent = L1TreeGetCenCode(targetComponent,pred,currentCenCode,rootID,nextID);
        rootID = optConnectionStrategy.Connection(end);
        treeGraph = graph(lastNeighborMatrix(:,1),lastNeighborMatrix(:,2),lastNeighborMatrix(:,3));
        [~,pred] = minspantree(treeGraph,'Method','dense','Type','tree','Root',rootID);
        nextID = find(pred == rootID);
        lastComponent = [];
        currentCenCode = 'T';
        lastComponent = L1TreeGetCenCode(lastComponent,pred,currentCenCode,rootID,nextID);
        lastComponentID = optConnectionStrategy.CurrentComponentID;
        for i = 1:1:size(lastComponent,2)
            for j = 1:1:length(lastComponent(i).SkeletonIDs)-1
                point_1 = allSkeletonPoints(lastComponent(i).SkeletonIDs(j),:);
                point_2 = allSkeletonPoints(lastComponent(i).SkeletonIDs(j+1),:);
                P1P2 = point_1 - point_2;
                P2P1 = point_2 - point_1;
                P1P0 = point_1 - tlsPoints;
                tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                tMatrix(tMatrix < 0) = 0;
                tMatrix(tMatrix > 1) = 1;
                pedalPoints = point_1 + tMatrix*P2P1;
                pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2));
                selectIDs = find(pedalDistances < tlsPointsToComponents(:,4));
                tlsPointsToComponents(selectIDs,1) = lastComponentID;
                tlsPointsToComponents(selectIDs,2) = i;
                tlsPointsToComponents(selectIDs,4) = pedalDistances(selectIDs);
            end
        end
        for i = 1:1:size(unfixedPoints,1)
            for j = 1:1:size(addSkeletonPoints,1)
                point_1 = unfixedPoints(i,:);
                point_2 = allSkeletonPoints(j,:);
                distance = sqrt(sum((point_1 - point_2).^2));
                if distance < parameters.minDistanceThreshold
                    pointID_1 = find((allSkeletonPoints(:,1) == point_1(1)) & (allSkeletonPoints(:,2) == point_1(2)) & (allSkeletonPoints(:,3) == point_1(3)));
                    pointID_2 = find((allSkeletonPoints(:,1) == point_2(1)) & (allSkeletonPoints(:,2) == point_2(2)) & (allSkeletonPoints(:,3) == point_2(3)));
                    unfixedNeighborMatrix = [unfixedNeighborMatrix;[pointID_1,pointID_2,distance]];
                end
            end
        end
        potentialConnections = graph(unfixedNeighborMatrix(:,1),unfixedNeighborMatrix(:,2),unfixedNeighborMatrix(:,3));
        delIDs = find([connectionStrategies(:).CurrentComponentID] == optConnectionStrategy.CurrentComponentID);
        connectionStrategies(delIDs) = [];
        delIDs = [];
        for i = 1:1:size(connectionStrategies,2)
            tempConnection = connectionStrategies(i).Connection;
            connectionParameters.CurrentComponentID = connectionStrategies(i).CurrentComponentID;
            selectID = find([Components(:).ID] == connectionStrategies(i).CurrentComponentID);
            connectionParameters.CurrentComponent = Components(selectID).Component;
            connectionParameters.StartDirection = connectionStrategies(i).StartDirection;
            connectionParameters.LastAddTLSPointIDs = connectionStrategies(i).LastAddTLSPointIDs;
            connectionParameters.LastSearchRange = connectionStrategies(i).LastSearchRange;
            [isCredible,connectionStrategy] = L1TreeIsConnectionCredible(tempConnection,targetComponent,allSkeletonPoints,tlsPoints,ignoreDirections,tlsPointsToComponents,connectionParameters,parameters);
            if isCredible
                connectionStrategies(i) = connectionStrategy;
            else
                delIDs = [delIDs,i];
            end
        end
        connectionStrategies(delIDs) = [];
        delID = find([Components(:).ID] == optConnectionStrategy.CurrentComponentID);
        Components(delID) = [];
        waitbar(0,globalBar,num2str(size(Components,2)));
    end
    % Step 4.5: Update iteration conditions
    if (size(Components,2) > oriComponentNum/10) && (parameters.maxDistanceThreshold < 2)
        selectIDs = find(~ismember(tlsPointsToComponents(:,1),[Components(:).ID]));
        tlsPointsToComponents(selectIDs,1) = 1;
        for i = 1:1:size(targetComponent,2)
            for j = 1:1:length(targetComponent(i).SkeletonIDs)-1
                point_1 = allSkeletonPoints(targetComponent(i).SkeletonIDs(j),:);
                point_2 = allSkeletonPoints(targetComponent(i).SkeletonIDs(j+1),:);
                P1P2 = point_1 - point_2;
                P2P1 = point_2 - point_1;
                P1P0 = point_1 - tlsPoints;
                tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                tMatrix(tMatrix < 0) = 0;
                tMatrix(tMatrix > 1) = 1;
                pedalPoints = point_1 + tMatrix*P2P1;
                pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
                selectIDs = find((tlsPointsToComponents(:,1) == 1) & (pedalDistances < tlsPointsToComponents(:,4)));
                tlsPointsToComponents(selectIDs,2) = i;
                tlsPointsToComponents(selectIDs,4) = pedalDistances(selectIDs);
            end
        end
        lastComponent = targetComponent;
        lastComponentID = 1;
        parameters.maxDistanceThreshold = parameters.maxDistanceThreshold + 0.2;
    else
        break;
    end
end
% Step 4.6: Build 3D tree model
for i = 1:1:size(targetComponent,2)
    delIDs = [];
    for j = 2:1:length(targetComponent(i).SkeletonIDs)-1
        if ismember(targetComponent(i).SkeletonIDs(j),allInsertPointIDs)
            delIDs = [delIDs,j];
        end
    end
    targetComponent(i).SkeletonIDs(delIDs) = [];
end
rootPoint = allSkeletonPoints(targetComponent(1).SkeletonIDs(1),:);
basePoint = parameters.basePoint;
currentDirection = allSkeletonPoints(targetComponent(1).SkeletonIDs(end),:) - allSkeletonPoints(targetComponent(1).SkeletonIDs(1),:);
currentDirection = currentDirection/sqrt(sum(currentDirection.^2));
newDirection = rootPoint - basePoint;
newDirection = newDirection/sqrt(sum(newDirection.^2));
angle = acos(currentDirection*newDirection')*180/pi;
if angle < 90
    allSkeletonPoints = [allSkeletonPoints;basePoint];
    targetComponent(1).SkeletonIDs = [size(allSkeletonPoints,1),targetComponent(1).SkeletonIDs];
end
treeModel = L1TreeUpdate(targetComponent,allSkeletonPoints,0);

%% Step 5: Remove unreasonable branches and update the 3D tree model
[treeModel,removeBranches] = L1TreeRemoveBranches(treeModel,tlsPoints,parameters);
if ~isempty(removeBranches)
    tempPoints = [endPoints(:).P];
    endPoints = reshape(tempPoints,3,size(endPoints,2))';
    delBranchIDs = [];
    for i = 1:1:size(removeBranches,2)
        tempPoint = removeBranches(i).Curve(1,:);
        for j = 1:1:size(treeModel,2)
            if ~treeModel(j).IsTerminal
                continue;
            end
            selectID = find((treeModel(j).Curve(:,1) == tempPoint(1)) & (treeModel(j).Curve(:,2) == tempPoint(2)) & (treeModel(j).Curve(:,3) == tempPoint(3)));
            if isempty(selectID) || (selectID < 3) || (selectID > size(treeModel(j).Curve,1)-2)
                continue;
            end
            consideredPoint_1 = treeModel(j).Curve(selectID-2,:);
            consideredPoint_2 = treeModel(j).Curve(selectID-1,:);
            consideredPoint_3 = treeModel(j).Curve(selectID+1,:);
            direction_1 = (consideredPoint_2 - consideredPoint_1)/sqrt(sum((consideredPoint_2 - consideredPoint_1).^2));
            direction_2 = (tempPoint - consideredPoint_2)/sqrt(sum((tempPoint - consideredPoint_2).^2));
            direction_3 = (consideredPoint_3 - tempPoint)/sqrt(sum((consideredPoint_3 - tempPoint).^2));
            angle_1 = acos(direction_1*direction_2')*180/pi;
            angle_2 = acos(direction_2*direction_3')*180/pi;
            isExist = find((endPoints(:,1) == consideredPoint_2(1)) & (endPoints(:,2) == consideredPoint_2(2)) & (endPoints(:,3) == consideredPoint_2(3)));
            if ~isempty(isExist) && (angle_1 + angle_2 > 150)
                delBranchIDs = [delBranchIDs;[j,selectID-1]];
            end
        end
    end
    if ~isempty(delBranchIDs)
        newDelBranchIDs = [];
        branchIDs = unique(delBranchIDs(:,1));
        for i = 1:1:length(branchIDs)
            selectIDs = find(delBranchIDs(:,1) == branchIDs(i));
            newDelBranchIDs = [newDelBranchIDs;[branchIDs(i),min(delBranchIDs(selectIDs,2))]];
        end
        delBranchIDs = newDelBranchIDs;
        newDelBranchIDs = [];
        for i = 1:1:size(delBranchIDs,1)
            tempRemoveBranch = treeModel(delBranchIDs(i,1));
            tempRemoveBranch.SkeletonIDs = tempRemoveBranch.SkeletonIDs(delBranchIDs(i,2):end);
            tempRemoveBranch.Curve = tempRemoveBranch.Curve(delBranchIDs(i,2):end,:);
            removeBranches = [removeBranches,tempRemoveBranch];
            treeModel(delBranchIDs(i,1)).SkeletonIDs = treeModel(delBranchIDs(i,1)).SkeletonIDs(1:delBranchIDs(i,2));
            treeModel(delBranchIDs(i,1)).Curve = treeModel(delBranchIDs(i,1)).Curve(1:delBranchIDs(i,2),:);
            for j = 1:1:size(treeModel,2)
                if contains(treeModel(j).CenCode,treeModel(delBranchIDs(i,1)).CenCode) && (treeModel(j).CenOrder > treeModel(delBranchIDs(i,1)).CenOrder)
                    newDelBranchIDs = [newDelBranchIDs,j];
                end
            end
        end
        removeBranches = [removeBranches,treeModel(newDelBranchIDs)];
        treeModel(newDelBranchIDs) = [];
    end
end

%% Step 6: Identify TLS points associated with each branch
% Step 6.1: Identify the removed branches
newRemoveBranches = [];
count = 1;
for i = 1:1:size(removeBranches,2)
    newRemoveBranches(count).Curve = removeBranches(i).Curve;
    count = count + 1;
end
for i = 1:1:size(Components,2)
    tempComponent = Components(i).Component;
    for j = 1:1:size(tempComponent,2)
        newRemoveBranches(count).Curve = allSkeletonPoints(tempComponent(j).SkeletonIDs,:);
        count = count + 1;
    end
end
% Step 6.2: Calculate the distances between TLS points and removed branches
tlsPointsToRemovedBranches = ones(size(tlsPoints,1),1)*9999;
for i = 1:1:size(newRemoveBranches,2)
    for j = 1:1:size(newRemoveBranches(i).Curve,1)-1
        point_1 = newRemoveBranches(i).Curve(j,:);
        point_2 = newRemoveBranches(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - tlsPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
        selectIDs = find(pedalDistances < tlsPointsToRemovedBranches);
        tlsPointsToRemovedBranches(selectIDs) = pedalDistances(selectIDs);
    end
end
% Step 6.3: Calculate the distances between unfixed skeleton points and branches, and the distances between TLS points and branches
unfixedToBranches = ones(size(unfixedPoints,1),2)*9999;
tlsPointsToBranches = ones(size(tlsPoints,1),3)*9999;
for i = 1:1:size(treeModel,2)
    for j = 1:1:size(treeModel(i).Curve,1)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        % Calculate the distances between unfixed skeleton points and branches
        P1P0 = point_1 - unfixedPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - unfixedPoints).^2,2));
        selectIDs = find(pedalDistances < unfixedToBranches(:,2));
        if length(find([treeModel(:).SkeletonIDs] == treeModel(i).SkeletonIDs(j))) == 1
            if j == 1
                selectIDs_1 = find((tMatrix == 0) & (pedalDistances < unfixedToBranches(:,2)));
                selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < unfixedToBranches(:,2)));
                unfixedToBranches(selectIDs_1,1) = 1;
                unfixedToBranches(selectIDs_2,1) = 0;
            elseif j == size(treeModel(i).Curve,1)-1
                selectIDs_1 = find((tMatrix == 1) & (pedalDistances < unfixedToBranches(:,2)));
                selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < unfixedToBranches(:,2)));
                unfixedToBranches(selectIDs_1,1) = 1;
                unfixedToBranches(selectIDs_2,1) = 0;
            else
                unfixedToBranches(selectIDs,1) = 0;
            end
        else
            unfixedToBranches(selectIDs,1) = 0;
        end
        unfixedToBranches(selectIDs,2) = pedalDistances(selectIDs);
        % Calculate the distances between TLS points and branches
        P1P0 = point_1 - tlsPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
        selectIDs = find(pedalDistances < tlsPointsToBranches(:,3));
        tlsPointsToBranches(selectIDs,1) = i;
        if length(find([treeModel(:).SkeletonIDs] == treeModel(i).SkeletonIDs(j))) == 1
            if j == 1
                selectIDs_1 = find((tMatrix == 0) & (pedalDistances < tlsPointsToBranches(:,3)));
                selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < tlsPointsToBranches(:,3)));
                tlsPointsToBranches(selectIDs_1,2) = 1;
                tlsPointsToBranches(selectIDs_2,2) = 0;
            elseif j == size(treeModel(i).Curve,1)-1
                selectIDs_1 = find((tMatrix == 1) & (pedalDistances < tlsPointsToBranches(:,3)));
                selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < tlsPointsToBranches(:,3)));
                tlsPointsToBranches(selectIDs_1,2) = 1;
                tlsPointsToBranches(selectIDs_2,2) = 0;
            else
                tlsPointsToBranches(selectIDs,2) = 0;
            end
        else
            tlsPointsToBranches(selectIDs,2) = 0;
        end
        tlsPointsToBranches(selectIDs,3) = pedalDistances(selectIDs);
    end
end
% Step 6.4: Identify associated and unassociated TLS points
tlsPoints = [tlsPoints,ones(size(tlsPoints,1),1)*9999];
for i = 1:1:size(treeModel,2)
    selectIDs = find((tlsPointsToBranches(:,1) == i) & (tlsPointsToBranches(:,2) == 0) & (tlsPointsToBranches(:,3) < tlsPointsToRemovedBranches));
    tempSearchRange = quantile(tlsPointsToBranches(selectIDs,3),0.9);
    delIDs_1 = find(unfixedToBranches(:,2) <= parameters.resolution);
    delIDs_2 = find((unfixedToBranches(:,1) == 0) & (unfixedToBranches(:,2) <= tempSearchRange));
    delIDs = unique([delIDs_1;delIDs_2]);
    currentUnfixedPoints = unfixedPoints;
    currentUnfixedPoints(delIDs,:) = [];
    if isempty(currentUnfixedPoints)
        selectIDs = find((tlsPointsToBranches(:,1) == i) & (tlsPointsToBranches(:,2) == 0) & (tlsPointsToBranches(:,3) < tlsPointsToRemovedBranches));
    else
        tlsPointsToUnfixedPoints = pdist2(currentUnfixedPoints,tlsPoints(:,1:3),'euclidean','Smallest',1);
        selectIDs = find((tlsPointsToBranches(:,1) == i) & (tlsPointsToBranches(:,2) == 0) & (tlsPointsToBranches(:,3) < tlsPointsToUnfixedPoints') & (tlsPointsToBranches(:,3) < tlsPointsToRemovedBranches));
    end
    tlsPoints(selectIDs,4) = i;
end
selectIDs = find(tlsPoints(:,4) == 9999);
tlsPoints(selectIDs,4) = 0;
baseBranchIDs = find([treeModel(:).CenOrder] <= 3);
otherBranchIDs = find([treeModel(:).CenOrder] > 3);
selectIDs = find(ismember(tlsPointsToBranches(:,1),baseBranchIDs) & (tlsPointsToBranches(:,3) > parameters.baseRadius*2));
tlsPoints(selectIDs,4) = 0;
selectIDs = find(ismember(tlsPointsToBranches(:,1),otherBranchIDs) & (tlsPointsToBranches(:,3) > parameters.baseRadius));
tlsPoints(selectIDs,4) = 0;

end