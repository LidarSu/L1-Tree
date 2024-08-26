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
% L1TreeSmoothBranches.m    The function for obtaining branch radii and lengths under different coding strategies
%
% Version 1.0
% Latest update     15 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% fileName    Filename for point cloud data without extensions
% treeModel    The input 3D tree model
% tlsPoints    The coordinates and attributes of TLS points
% globalBar    The global waiting bar
%
% OUTPUTS:
% treeModelCentrifugal    The 3D tree model under the centrifugal coding strategy
% treeModelGaaliche    The 3D tree model under Gaaliche's coding strategy
% ------------------------------------------------------------------------------

function [treeModelCentrifugal,treeModelGaaliche] = L1TreeConstructModel(fileName,treeModel,tlsPoints,globalBar)

waitbar(5/5,globalBar,'Construct the 3D tree model');

%% Step 1: Remove anomalies in each branch
for i = 1:1:size(treeModel,2)
    if size(treeModel(i).Curve,1) < 3
        continue;
    end
    newCurve = treeModel(i).Curve;
    while true
        angles = [];
        for j = 2:1:size(newCurve,1)-1
            point_1 = newCurve(j-1,:);
            point_2 = newCurve(j,:);
            point_3 = newCurve(j+1,:);
            direction_1 = point_2 - point_1;
            direction_1 = direction_1/sqrt(sum(direction_1.^2));
            direction_2 = point_3 - point_2;
            direction_2 = direction_2/sqrt(sum(direction_2.^2));
            angles = [angles;acos(direction_1*direction_2')*180/pi];
        end
        if all(angles < 90)
            break;
        end
        [~,maxID] = max(angles);
        newCurve(maxID+1,:) = [];
    end
    treeModel(i).Curve = newCurve;
end
skeletonPoints = [];
for i = 1:1:size(treeModel,2)
    for j = 1:1:size(treeModel(i).Curve,1)
        tempPoint = treeModel(i).Curve(j,:);
        if isempty(skeletonPoints)
            skeletonPoints = [skeletonPoints;tempPoint];
            continue;
        end
        selectID = find((skeletonPoints(:,1) == tempPoint(1)) & (skeletonPoints(:,2) == tempPoint(2)) & (skeletonPoints(:,3) == tempPoint(3)));
        if isempty(selectID)
            skeletonPoints = [skeletonPoints;tempPoint];
        end
    end
end
for i = 1:1:size(treeModel,2)
    newSkeletonIDs = [];
    for j = 1:1:size(treeModel(i).Curve,1)
        tempPoint = treeModel(i).Curve(j,:);
        selectID = find((skeletonPoints(:,1) == tempPoint(1)) & (skeletonPoints(:,2) == tempPoint(2)) & (skeletonPoints(:,3) == tempPoint(3)));
        newSkeletonIDs = [newSkeletonIDs,selectID];
    end
    treeModel(i).SkeletonIDs = newSkeletonIDs;
end

%% Step 2: Update the associated branch ID of ID of each TLS point
lengths = ones(size(treeModel,2),1)*9999;
for i = 1:1:size(treeModel,2)
    totalLength = 0;
    for j = 1:1:size(treeModel(i).Curve,1)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        segmentLength = sqrt(sum((point_1 - point_2).^2));
        totalLength = totalLength + segmentLength;
    end
    lengths(i) = totalLength;
end
tlsPointsToBranches = ones(size(tlsPoints,1),3)*9999;
for i = 1:1:size(treeModel,2)
    for j = 1:1:size(treeModel(i).Curve,1)-1
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
maxCenOrder = max([treeModel(:).CenOrder]);
associatedTLSPointIDs = [];
for i = 1:1:maxCenOrder
    for j = 1:1:size(treeModel,2)
        if treeModel(j).CenOrder ~= i
            continue;
        end
        if i == 1
            selectIDs = find((tlsPointsToBranches(:,1) == j) & (tlsPointsToBranches(:,2) == 0));
            searchRange = quantile(tlsPointsToBranches(selectIDs,3),0.9)*2;
        else
            selectIDs = find((tlsPointsToBranches(:,1) == j) & (tlsPointsToBranches(:,2) == 0));
            searchRange = quantile(tlsPointsToBranches(selectIDs,3),0.9)*2;
            for k = 1:1:size(treeModel,2)
                if contains(treeModel(j).CenCode,treeModel(k).CenCode) && (k ~= j) && (searchRanges(k) < searchRange) && (lengths(k) > 0.2)
                    searchRange = searchRanges(k);
                end
            end
        end
        searchRanges(j) = searchRange;
        for k = 1:1:size(treeModel(j).Curve,1)-1
            point_1 = treeModel(j).Curve(k,:);
            point_2 = treeModel(j).Curve(k+1,:);
            P1P2 = point_1 - point_2;
            P2P1 = point_2 - point_1;
            P1P0 = point_1 - tlsPoints(:,1:3);
            tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
            tMatrix(tMatrix < 0) = 0;
            tMatrix(tMatrix > 1) = 1;
            pedalPoints = point_1 + tMatrix*P2P1;
            pedalDistances = sqrt(sum((pedalPoints - tlsPoints(:,1:3)).^2,2));
            selectIDs = find((tlsPointsToBranches(:,2) == 0) & (pedalDistances <= searchRange));
            associatedTLSPointIDs = [associatedTLSPointIDs;selectIDs];
        end
    end
end
associatedTLSPointIDs = sort(unique(associatedTLSPointIDs));
selectIDs_1 = find(tlsPoints(:,4) == 0);
selectIDs_2 = setdiff(1:size(tlsPoints,1),associatedTLSPointIDs);
tlsPoints(selectIDs_1,4) = 0;
tlsPoints(selectIDs_2,4) = 0;

%% Step 3: Calculate branch radii under the centrifugal coding strategy
for i = 1:1:maxCenOrder
    for j = 1:1:size(treeModel,2)
        if treeModel(j).CenOrder ~= i
            continue;
        end
        selectIDs = find(tlsPoints(:,4) == j);
        if isempty(selectIDs)
            treeModel(j).CenRadius = nan;
            continue;
        end
        currentTLSPoints = tlsPoints(selectIDs,1:3);
        tlsPointsToBranch = ones(size(currentTLSPoints,1),2)*9999;
        for k = 1:1:size(treeModel(j).Curve,1)-1
            point_1 = treeModel(j).Curve(k,:);
            point_2 = treeModel(j).Curve(k+1,:);
            P1P2 = point_1 - point_2;
            P2P1 = point_2 - point_1;
            P1P0 = point_1 - currentTLSPoints;
            tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
            tMatrix(tMatrix < 0) = 0;
            tMatrix(tMatrix > 1) = 1;
            pedalPoints = point_1 + tMatrix*P2P1;
            pedalDistances = sqrt(sum((pedalPoints - currentTLSPoints).^2,2));
            selectIDs = find(pedalDistances < tlsPointsToBranch(:,2));
            if k == 1
                selectIDs_1 = find((tMatrix == 0) & (pedalDistances < tlsPointsToBranch(:,2)));
                selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < tlsPointsToBranch(:,2)));
                tlsPointsToBranch(selectIDs_1,1) = 1;
                tlsPointsToBranch(selectIDs_2,1) = 0;
            elseif k == size(treeModel(j).Curve,1)-1
                selectIDs_1 = find((tMatrix == 1) & (pedalDistances < tlsPointsToBranch(:,2)));
                selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < tlsPointsToBranch(:,2)));
                tlsPointsToBranch(selectIDs_1,1) = 1;
                tlsPointsToBranch(selectIDs_2,1) = 0;
            else
                tlsPointsToBranch(selectIDs,1) = 0;
            end
            tlsPointsToBranch(selectIDs,2) = pedalDistances(selectIDs);
        end
        if i == 1
            selectIDs = find(tlsPointsToBranch(:,1) == 0);
            treeModel(j).CenRadius = mean(tlsPointsToBranch(selectIDs,2));
            continue;
        end
        neighborBranchIDs = [];
        for k = 1:1:size(treeModel,2)
            if contains(treeModel(j).CenCode,treeModel(k).CenCode) && (treeModel(k).CenOrder == treeModel(j).CenOrder-1)
                neighborBranchIDs = [neighborBranchIDs,k];
            end
        end
        for k = 1:1:size(treeModel,2)
            if contains(treeModel(k).CenCode,treeModel(neighborBranchIDs(1)).CenCode) && (treeModel(k).CenOrder == i) && (k ~= j)
                neighborBranchIDs = [neighborBranchIDs,k];
            end
        end
        parentRadius = treeModel(neighborBranchIDs(1)).CenRadius*2;
        delIDs = [];
        for k = 1:1:length(neighborBranchIDs)
            tlsPointsToNeighborBranch = ones(size(currentTLSPoints,1),2)*9999;
            for l = 1:1:size(treeModel(neighborBranchIDs(k)).Curve,1)-1
                point_1 = treeModel(neighborBranchIDs(k)).Curve(l,:);
                point_2 = treeModel(neighborBranchIDs(k)).Curve(l+1,:);
                P1P2 = point_1 - point_2;
                P2P1 = point_2 - point_1;
                P1P0 = point_1 - currentTLSPoints;
                tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                tMatrix(tMatrix < 0) = 0;
                tMatrix(tMatrix > 1) = 1;
                pedalPoints = point_1 + tMatrix*P2P1;
                pedalDistances = sqrt(sum((pedalPoints - currentTLSPoints).^2,2));
                selectIDs = find(pedalDistances < tlsPointsToNeighborBranch(:,2));
                if l == 1
                    selectIDs_1 = find((tMatrix == 0) & (pedalDistances < tlsPointsToNeighborBranch(:,2)));
                    selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < tlsPointsToNeighborBranch(:,2)));
                    tlsPointsToNeighborBranch(selectIDs_1,1) = 1;
                    tlsPointsToNeighborBranch(selectIDs_2,1) = 0;
                elseif l == size(treeModel(neighborBranchIDs(k)).Curve,1)-1
                    selectIDs_1 = find((tMatrix == 1) & (pedalDistances < tlsPointsToNeighborBranch(:,2)));
                    selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < tlsPointsToNeighborBranch(:,2)));
                    tlsPointsToNeighborBranch(selectIDs_1,1) = 1;
                    tlsPointsToNeighborBranch(selectIDs_2,1) = 0;
                else
                    tlsPointsToNeighborBranch(selectIDs,1) = 0;
                end
                tlsPointsToNeighborBranch(selectIDs,2) = pedalDistances(selectIDs);
            end
            delIDs = [delIDs;find(tlsPointsToNeighborBranch(:,2) < parentRadius)];
        end
        delIDs = unique(delIDs);
        tlsPointsToBranch(delIDs,:) = [];
        selectIDs = find(tlsPointsToBranch(:,1) == 0);
        treeModel(j).CenRadius = mean(tlsPointsToBranch(selectIDs,2));
    end
end
maxCenOrder = max([treeModel(:).CenOrder]);
incredibleBranchIDs = [];
for i = 2:1:maxCenOrder
    for j = 1:1:size(treeModel,2)
        if treeModel(j).CenOrder ~= i
            continue;
        end
        parentBranchID = [];
        childBranchIDs = [];
        for k = 1:1:size(treeModel,2)
            if contains(treeModel(j).CenCode,treeModel(k).CenCode) && (treeModel(k).CenOrder == i-1)
                parentBranchID = k;
            end
            if contains(treeModel(k).CenCode,treeModel(j).CenCode) && (treeModel(k).CenOrder == i + 1) && ~isnan(treeModel(k).CenRadius)
                childBranchIDs = [childBranchIDs,k];
            end
        end
        if isempty(childBranchIDs) && (isnan(treeModel(j).CenRadius) || (treeModel(j).CenRadius > treeModel(parentBranchID).CenRadius))
            treeModel(j).CenRadius = treeModel(parentBranchID).CenRadius;
            incredibleBranchIDs = [incredibleBranchIDs,j];
            continue;
        end
        if ~isempty(childBranchIDs) && (isnan(treeModel(j).CenRadius) || (treeModel(j).CenRadius > treeModel(parentBranchID).CenRadius))
            parentRadius = treeModel(parentBranchID).CenRadius;
            childRadii = [treeModel(childBranchIDs).CenRadius];
            if max(childRadii) > parentRadius
                treeModel(j).CenRadius = parentRadius;
                incredibleBranchIDs = [incredibleBranchIDs,j];
            else
                treeModel(j).CenRadius = (max(childRadii) + parentRadius)/2;
            end
            continue;
        end
    end
end
incredibleBranchIDs = unique(incredibleBranchIDs);
for i = 1:1:length(incredibleBranchIDs)
    treeModel(incredibleBranchIDs(i)).IsCenCredible = false;
    treeModel(incredibleBranchIDs(i)).IsGaCredible = false;
end

%% Step 4: Calculate branch lengths under the centrifugal coding strategy
for i = 1:1:size(treeModel,2)
    totalLength = 0;
    for j = 1:1:size(treeModel(i).Curve,1)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        segmentLength = sqrt(sum((point_1 - point_2).^2));
        totalLength = totalLength + segmentLength;
    end
    treeModel(i).CenLength = totalLength;
end

%% Step 5: Output the 3D tree model under the centrifugal coding strategy
for i = 1:1:size(treeModel,2)
    treeModelCentrifugal(i).Curve = treeModel(i).Curve;
    treeModelCentrifugal(i).CenCode = treeModel(i).CenCode;
    treeModelCentrifugal(i).CenOrder = treeModel(i).CenOrder;
    treeModelCentrifugal(i).CenRadius = treeModel(i).CenRadius;
    treeModelCentrifugal(i).CenLength = treeModel(i).CenLength;
    treeModelCentrifugal(i).IsCenCredible = treeModel(i).IsCenCredible;
end
outputFileName = ['./Results/',fileName,'_Centrifugal.mat'];
save(outputFileName,'treeModelCentrifugal');

%% Step 6: Obtain transit codes
maxCenOrder = max([treeModel(:).CenOrder]);
for i = 1:1:maxCenOrder
    for j = 1:1:size(treeModel,2)
        if (treeModel(j).CenOrder ~= i) || ~isempty(treeModel(j).GaCode)
            continue;
        end
        if i == 1
            treeModel(j).GaCode = 'T';
            break;
        end
        for k = 1:1:size(treeModel,2)
            if contains(treeModel(j).CenCode,treeModel(k).CenCode) && (treeModel(k).CenOrder == treeModel(j).CenOrder-1)
                parentBranchID = k;
                break;
            end
        end
        brotherBranchIDs = [];
        for k = 1:1:size(treeModel,2)
            if contains(treeModel(k).CenCode,treeModel(parentBranchID).CenCode) && (treeModel(k).CenOrder == treeModel(parentBranchID).CenOrder+1)
                brotherBranchIDs = [brotherBranchIDs,k];
            end
        end
        summaryTable = [];
        parentDirection = treeModel(parentBranchID).Curve(end,:) - treeModel(parentBranchID).Curve(1,:);
        parentDirection = parentDirection/sqrt(sum(parentDirection.^2));
        for k = 1:1:length(brotherBranchIDs)
            tempDirection = treeModel(brotherBranchIDs(k)).Curve(end,:) - treeModel(brotherBranchIDs(k)).Curve(1,:);
            tempDirection = tempDirection/sqrt(sum(tempDirection.^2));
            angle = acos(tempDirection*parentDirection')*180/pi;
            tempRadius = treeModel(brotherBranchIDs(k)).CenRadius;
            summaryTable = [summaryTable;[brotherBranchIDs(k),tempRadius,angle]];
        end
        % Identify the continuation of the current parent branch
        if isempty(find(isnan(summaryTable(:,2))))
            summaryTable = sortrows(summaryTable,3,'ascend');
            continueBranchID = summaryTable(1,1);
        else
            summaryTable = sortrows(summaryTable,2,'descend');
            minAngle = min(summaryTable(:,3));
            for k = 1:1:size(summaryTable,1)
                if summaryTable(k,3)/minAngle <= 1.1
                    continueBranchID = summaryTable(k,1);
                    break;
                end
            end
        end
        count = 1;
        for k = 1:1:size(summaryTable,1)
            if summaryTable(k,1) == continueBranchID
                treeModel(summaryTable(k,1)).GaCode = [treeModel(parentBranchID).GaCode,'_C'];
            else
                treeModel(summaryTable(k,1)).GaCode = [treeModel(parentBranchID).GaCode,'_B',num2str(count)];
                count = count + 1;
            end
        end
    end
end

%% Step 7: Obtain branch codes, orders, radii, and lengths under Gaaliche's coding strategy 
% based on transit codes and the atrributes under the centrifugal coding strategy
treeModelGaaliche = [];
flag = false(1,size(treeModel,2));
while any(~flag)
    selectID = find(~flag,1);
    rootGaCode = treeModel(selectID).GaCode;
    while strcmp(rootGaCode(end),'C')
        rootGaCode = rootGaCode(1:(end-2));
    end
    for i = 1:1:size(treeModel,2)
        if strcmp(treeModel(i).GaCode,rootGaCode)
            rootBranchID = i;
            break;
        end
    end
    associatedBranchIDs = [];
    currentGaCode = rootGaCode;
    while true
        isExist = false;
        for i = 1:1:size(treeModel,2)
            if strcmp(treeModel(i).GaCode,currentGaCode)
                associatedBranchIDs = [associatedBranchIDs,i];
                isExist = true;
            end
        end
        if isExist
            currentGaCode = [currentGaCode,'_C'];
        else
            break;
        end
    end
    % Obtain the attributes of a branch under Gaaliche's coding strategy
    newBranch.Curve = [];
    for i = 1:1:length(associatedBranchIDs)
        if i == 1
            newBranch.Curve = treeModel(associatedBranchIDs(1)).Curve;
        else
            newBranch.Curve = [newBranch.Curve;treeModel(associatedBranchIDs(i)).Curve(2:end,:)];
        end
    end
    newBranch.GaCode = rootGaCode;
    TNum = length(strfind(newBranch.GaCode,'T'));
    BNum = length(strfind(newBranch.GaCode,'B'));
    newBranch.GaOrder = TNum + BNum;
    newBranch.GaRadius = treeModel(rootBranchID).CenRadius;
    newBranch.GaLength = 0;
    for i = 1:1:size(newBranch.Curve,1)-1
        point_1 = newBranch.Curve(i,:);
        point_2 = newBranch.Curve(i+1,:);
        segmentLength = sqrt(sum((point_1 - point_2).^2));
        newBranch.GaLength = newBranch.GaLength + segmentLength;
    end
    newBranch.IsGaCredible = treeModel(rootBranchID).IsGaCredible;
    treeModelGaaliche = [treeModelGaaliche,newBranch];
    flag(associatedBranchIDs) = true;
end

%% Step 8: Output the 3D tree model under Gaaliche's coding strategy
outputFileName = ['./Results/',fileName,'_Gaaliche.mat'];
save(outputFileName,'treeModelGaaliche');

end