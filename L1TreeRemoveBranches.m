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
% L1TreeRemoveBranches.m    The function for removing unreasonable branches
%
% Version 1.0
% Latest update     04 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% treeModel    The input 3D tree model
% tlsPoints    The coordinates of TLS points
% parameters    The parameter set
%
% OUTPUTS:
% treeModel    The update 3D tree model
% removeBranches    The branches need to be removed
% ------------------------------------------------------------------------------

function [treeModel,removeBranches] = L1TreeRemoveBranches(treeModel,tlsPoints,parameters)

%% Step 1: Obtain the coordinates of skeleton points
skeletonPoints = [];
for i = 1:1:size(treeModel,2)
    skeletonPoints = [skeletonPoints;[treeModel(i).SkeletonIDs',treeModel(i).Curve]];
end
skeletonPoints = unique(skeletonPoints,'rows','stable');
skeletonPoints = sortrows(skeletonPoints,1);
skeletonPoints = skeletonPoints(:,2:4);

%% Step 2: Apply length constraints
delBranchIDs = [];
for i = 1:1:size(treeModel,2)
    totalLength = 0;
    for j = 1:1:size(treeModel(i).Curve,1)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        totalLength = totalLength + sqrt(sum((point_1 - point_2).^2));
    end
    if (totalLength < parameters.resolution*(parameters.branchSizeThreshold-1)) && (length(treeModel(i).SkeletonIDs) < parameters.branchSizeThreshold)
        delBranchIDs = [delBranchIDs,i];
    end
end
for i = 1:1:length(delBranchIDs)
    if treeModel(delBranchIDs(i)).IsTerminal
        continue;
    end
    childBranchIDs = [];
    for j = 1:1:size(treeModel,2)
        if contains(treeModel(j).CenCode,treeModel(delBranchIDs(i)).CenCode) && (treeModel(j).CenOrder == treeModel(delBranchIDs(i)).CenOrder+1)
            childBranchIDs = [childBranchIDs,j];
        end
    end
    parentStartPointID = treeModel(delBranchIDs(i)).SkeletonIDs(1);
    parentStartPoint = treeModel(delBranchIDs(i)).Curve(1,:);
    for j = 1:1:length(childBranchIDs)
        treeModel(childBranchIDs(j)).SkeletonIDs = [parentStartPointID,treeModel(childBranchIDs(j)).SkeletonIDs(2:end)];
        treeModel(childBranchIDs(j)).Curve = [parentStartPoint;treeModel(childBranchIDs(j)).Curve(2:end,:)];
    end
end
treeModel(delBranchIDs) = [];
rootID = treeModel(1).SkeletonIDs(1);
neighborMatrix = [];
for i = 1:1:size(treeModel,2)
    for j = 1:1:length(treeModel(i).SkeletonIDs)-1
        pointID_1 = treeModel(i).SkeletonIDs(j);
        pointID_2 = treeModel(i).SkeletonIDs(j+1);
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        distance = sqrt(sum((point_1 - point_2).^2));
        neighborMatrix = [neighborMatrix;[pointID_1,pointID_2,distance]];
    end
end
[treeModel,skeletonPoints] = L1TreeUpdate(neighborMatrix,skeletonPoints,rootID);

%% Step 3: Apply angle constraints
delBranchIDs = [];
maxCenOrder = max([treeModel(:).CenOrder]);
for i = 1:1:maxCenOrder
    for j = 1:1:size(treeModel,2)
        if treeModel(j).CenOrder ~= i
            continue;
        end
        parentDirections = treeModel(j).Curve(2:end,:) - treeModel(j).Curve(1:end-1,:);
        parentDirections = parentDirections./sqrt(sum(parentDirections.^2,2));
        for k = 1:1:size(treeModel,2)
            if contains(treeModel(k).CenCode,treeModel(j).CenCode) && (treeModel(k).CenOrder == treeModel(j).CenOrder+1)
                childDirection = treeModel(k).Curve(end,:) - treeModel(k).Curve(1,:);
                childDirection = childDirection/sqrt(sum(childDirection.^2));
                angle = median(acos(parentDirections*childDirection')*180/pi);
                if angle > parameters.maxAngleThreshold
                    delBranchIDs = [delBranchIDs,k];
                end
            end
        end
    end
end
removeBranches = [];
if ~isempty(delBranchIDs)
    newDelBranchIDs = delBranchIDs;
    for i = 1:1:length(newDelBranchIDs)
        for j = 1:1:size(treeModel,2)
            if contains(treeModel(j).CenCode,treeModel(delBranchIDs(i)).CenCode) && (treeModel(j).CenOrder > treeModel(delBranchIDs(i)).CenOrder)
                newDelBranchIDs = [newDelBranchIDs,j];
            end
        end
    end
    newDelBranchIDs = unique(newDelBranchIDs);
    removeBranches = treeModel(newDelBranchIDs);
    treeModel(newDelBranchIDs) = [];
end

%% Step 4: Remove the branches that are within the influence of other branches
lengths = ones(size(treeModel,2),1)*9999;
for i = 1:1:size(treeModel,2)
    totalLength = 0;
    for j = 1:1:size(treeModel(i).Curve,1)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        totalLength = totalLength + sqrt(sum((point_1 - point_2).^2));
    end
    lengths(i) = totalLength;
end
[~,idx] = sort(lengths,'descend');
tlsPointsToBranches = ones(size(tlsPoints,1),3)*9999;
for i = 1:1:size(treeModel,2)
    for j = 1:1:size(treeModel(i).Curve,1)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
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
searchRanges = ones(size(treeModel,2),1)*9999;
maxCenOrder = max([treeModel(:).CenCode]);
for i = 1:1:maxCenOrder
    for j = 1:1:size(treeModel,2)
        if treeModel(j).CenOrder ~= i
            continue;
        end
        if i == 1
            selectIDs = find((tlsPointsToBranches(:,1) == j) & (tlsPointsToBranches(:,2) == 0));
            searchRanges(j) = median(tlsPointsToBranches(selectIDs,3))*2;
        else
            selectIDs = find((tlsPointsToBranches(:,1) == j) & (tlsPointsToBranches(:,2) == 0));
            tempSearchRange = quantile(tlsPointsToBranches(selectIDs,3),0.9)*2;
            for k = 1:1:size(treeModel,2)
                if contains(treeModel(j).CenCode,treeModel(k).CenCode) && (k ~= j) && (searchRanges(k) < tempSearchRange) && (lengths(k) > 0.2)
                    tempSearchRange = searchRanges(k);
                end
            end
            searchRanges(j) = tempSearchRange;
        end
    end
end
delBranchIDs = [];
for i = 1:1:length(idx)
    if ismember(idx(i),delBranchIDs)
        continue;
    end
    currentBranch = treeModel(idx(i));
    searchRange = searchRanges(idx(i));
    for j = 1:1:size(treeModel,2)
        if ismember(j,delBranchIDs) || (j == idx(i))
            continue;
        end
        consideredBranchPoints = treeModel(j).Curve;
        consideredToCurrentBranch = ones(size(consideredBranchPoints,1),2)*9999;
        for k = 1:1:size(currentBranch.Curve,1)-1
            point_1 = currentBranch.Curve(k,:);
            point_2 = currentBranch.Curve(k+1,:);
            P1P2 = point_1 - point_2;
            P2P1 = point_2 - point_1;
            P1P0 = point_1 - consideredBranchPoints;
            tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
            tMatrix(tMatrix < 1e-6) = 0;
            tMatrix(tMatrix > 1 - 1e-6) = 1;
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
        if (max(consideredToCurrentBranch(selectIDs,2)) < searchRange) && (length(selectIDs) > size(consideredToCurrentBranch,1)/2) && treeModel(j).IsTerminal
            delBranchIDs = [delBranchIDs,j];
        end
    end
end
treeModel(delBranchIDs) = [];

%% Step 5: Update the 3D tree model
rootID = treeModel(1).SkeletonIDs(1);
neighborMatrix = [];
for i = 1:1:size(treeModel,2)
    for j = 1:1:length(treeModel(i).SkeletonIDs)-1
        pointID_1 = treeModel(i).SkeletonIDs(j);
        pointID_2 = treeModel(i).SkeletonIDs(j+1);
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        distance = sqrt(sum((point_1 - point_2).^2));
        neighborMatrix = [neighborMatrix;[pointID_1,pointID_2,distance]];
    end
end
treeModel = L1TreeUpdate(neighborMatrix,skeletonPoints,rootID);

end