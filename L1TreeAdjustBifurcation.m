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
% L1TreeAdjustBifurcation.m    The function for adjusting bifurcation positions
%
% Version 1.0
% Latest update     06 Feb 2024
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

function [treeModel,tlsPoints] = L1TreeAdjustBifurcation(treeModel,tlsPoints,parameters,globalBar)

waitbar(2/5,globalBar,'Adjust bifurcation positions');

for i = 1:1:size(treeModel,2)
    treeModel(i).Curve = round(treeModel(i).Curve,6);
end

%% Step 1: Adjust bifurcation locations
newBifurcationPoints = [];    % Bifurcation points that have been optimized
while true
    % Step 1.1: Identify bifurcation points that have not been optimized
    bifurcationInfo = [];
    for i = 1:1:size(treeModel,2)
        bifurcationID = treeModel(i).SkeletonIDs(1);
        bifurcationPoint = treeModel(i).Curve(1,:);
        associatedBranchIDs = [];
        for j = 1:1:size(treeModel,2)
            if ismember(bifurcationID,treeModel(j).SkeletonIDs)
                associatedBranchIDs = [associatedBranchIDs,j];
            end
        end
        if length(associatedBranchIDs) >= 3
            bifurcationInfo = [bifurcationInfo;[bifurcationID,bifurcationPoint,min([treeModel(associatedBranchIDs).CenOrder])]];
        end
        bifurcationID = treeModel(i).SkeletonIDs(end);
        bifurcationPoint = treeModel(i).Curve(end,:);
        associatedBranchIDs = [];
        for j = 1:1:size(treeModel,2)
            if ismember(bifurcationID,treeModel(j).SkeletonIDs)
                associatedBranchIDs = [associatedBranchIDs,j];
            end
        end
        if length(associatedBranchIDs) >= 3
            bifurcationInfo = [bifurcationInfo;[bifurcationID,bifurcationPoint,min([treeModel(associatedBranchIDs).CenOrder])]];
        end
    end
    bifurcationInfo = unique(bifurcationInfo,'rows','stable');
    delIDs = [];
    for i = 1:1:size(bifurcationInfo,1)
        if isempty(newBifurcationPoints)
            break;
        end
        selectID = find((newBifurcationPoints(:,1) == bifurcationInfo(i,2)) & (newBifurcationPoints(:,2) == bifurcationInfo(i,3)) & (newBifurcationPoints(:,3) == bifurcationInfo(i,4)));
        if ~isempty(selectID)
            delIDs = [delIDs,i];
        end
    end
    bifurcationInfo(delIDs,:) = [];
    if isempty(bifurcationInfo)
        break;
    end
    % Step 1.2: Identify the branch point to be examined, the associated TLS points, and the ID of the assocaited parent branch
    selectID = find(bifurcationInfo(:,end) == min(bifurcationInfo(:,end)),1);
    bifurcationID = bifurcationInfo(selectID,1);
    bifurcationPoint = bifurcationInfo(selectID,2:4);
    associatedBranchIDs = [];
    for i = 1:1:size(treeModel,2)
        if ismember(bifurcationID,treeModel(i).SkeletonIDs)
            associatedBranchIDs = [associatedBranchIDs,i];
        end
    end
    associatedTLSPointIDs = find(ismember(tlsPoints(:,4),associatedBranchIDs));
    associatedTLSPoints = tlsPoints(associatedTLSPointIDs,1:3);
    minCenOrder = 9999;
    parentBranchID = 9999;
    for i = 1:1:length(associatedBranchIDs)
        if treeModel(associatedBranchIDs(i)).CenOrder < minCenOrder
            minCenOrder = treeModel(associatedBranchIDs(i)).CenOrder;
            parentBranchID = associatedBranchIDs(i);
        end
    end
    % Step 1.3: Identify all reasonable adjustment strategies
    adjustmentStrategies = [];
    for i = 1:1:length(associatedBranchIDs)
        targetBranch = treeModel(associatedBranchIDs(i));
        if targetBranch.SkeletonIDs(1) == bifurcationID
            startID = 2;
            endID = length(targetBranch.SkeletonIDs);
            stepLength = 1;
        end
        if targetBranch.SkeletonIDs(end) == bifurcationID
            startID = length(targetBranch.SkeletonIDs)-1;
            endID = 1;
            stepLength = -1;
        end
        for j = startID:stepLength:endID
            if j == startID
                lastPoint = bifurcationPoint;
            else
                lastPoint = targetBranch.Curve(j-stepLength,:);
            end
            potentialTargetPoint = targetBranch.Curve(j,:);
            segmentLength = sqrt(sum((potentialTargetPoint - lastPoint).^2));
            sampleNum = floor(segmentLength/parameters.resolution);
            for k = 1:1:(sampleNum+1)
                if k == sampleNum+1
                    isInterpolation = false;
                else
                    isInterpolation = true;
                end
                targetID = targetBranch.SkeletonIDs(j);
                targetPoint = lastPoint + (potentialTargetPoint - lastPoint)/(sampleNum+1)*k;
                targetPoint = round(targetPoint,6);
                if (targetBranch.Curve(endID,1) == targetPoint(1)) && (targetBranch.Curve(endID,2) == targetPoint(2)) && (targetBranch.Curve(endID,3) == targetPoint(3))
                    isTerminal = true;
                else
                    isTerminal = false;
                end
                % Affected branch segments (), unaffected branch segments, and newly add branch segments
                branchSegments_1 = [];
                for l = 1:1:length(associatedBranchIDs)
                    tempBranch = treeModel(associatedBranchIDs(l));
                    selectID = find(tempBranch.SkeletonIDs == bifurcationID);
                    if (associatedBranchIDs(l) == associatedBranchIDs(i)) && (selectID == 1)
                        for m = 1:1:(j-2)
                            point_1 = tempBranch.Curve(m,:);
                            point_2 = tempBranch.Curve(m+1,:);
                            branchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                        if isInterpolation
                            point_1 = tempBranch.Curve(j-1,:);
                            point_2 = targetPoint;
                            branchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                        else
                            point_1 = tempBranch.Curve(j-1,:);
                            point_2 = tempBranch.Curve(j,:);
                            ranchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                    end
                    if (associatedBranchIDs(l) == associatedBranchIDs(i)) && (selectID ~= 1)
                        for m = length(tempBranch.SkeletonIDs):-1:(j+2)
                            point_1 = tempBranch.Curve(m,:);
                            point_2 = tempBranch.Curve(m-1,:);
                            branchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                        if isInterpolation
                            point_1 = tempBranch.Curve(j+1,:);
                            point_2 = targetPoint;
                            branchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                        else
                            point_1 = tempBranch.Curve(j+1,:);
                            point_2 = tempBranch.Curve(j,:);
                            branchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                    end
                    if (associatedBranchIDs(l) ~= associatedBranchIDs(i)) && (selectID == 1)
                        point_1 = tempBranch.Curve(1,:);
                        point_2 = tempBranch.Curve(2,:);
                        branchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                    end
                    if (associatedBranchIDs(l) ~= associatedBranchIDs(i)) && (selectID ~= 1)
                        point_1 = tempBranch.Curve(end,:);
                        point_2 = tempBranch.Curve(end-1,:);
                        branchSegments_1 = [branchSegments_1;[point_1,point_2,associatedBranchIDs(l)]];
                    end
                end
                branchSegments_2 = [];
                for l = 1:1:length(associatedBranchIDs)
                    tempBranch = treeModel(associatedBranchIDs(l));
                    selectID = find(tempBranch.SkeletonIDs == bifurcationID);
                    if (associatedBranchIDs(l) == associatedBranchIDs(i)) && (selectID == 1)
                        if isInterpolation
                            point_1 = targetPoint;
                            point_2 = tempBranch.Curve(j,:);
                            branchSegments_2 = [branchSegments_2;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                        for m = j:1:length(tempBranch.SkeletonIDs)-1
                            point_1 = tempBranch.Curve(m,:);
                            point_2 = tempBranch.Curve(m+1,:);
                            branchSegments_2 = [branchSegments_2;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                    end
                    if (associatedBranchIDs(l) == associatedBranchIDs(i)) && (selectID ~= 1)
                        if isInterpolation
                            point_1 = targetPoint;
                            point_2 = tempBranch.Curve(j,:);
                            branchSegments_2 = [branchSegments_2;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                        for m = j:-1:2
                            point_1 = tempBranch.Curve(m,:);
                            point_2 = tempBranch.Curve(m-1,:);
                            branchSegments_2 = [branchSegments_2;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                    end
                    if (associatedBranchIDs(l) ~= associatedBranchIDs(i)) && (selectID == 1)
                        for m = 2:1:length(tempBranch.SkeletonIDs)-1
                            point_1 = tempBranch.Curve(m,:);
                            point_2 = tempBranch.Curve(m+1,:);
                            branchSegments_2 = [branchSegments_2;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                    end
                    if (associatedBranchIDs(l) ~= associatedBranchIDs(i)) && (selectID ~= 1)
                        for m = (length(tempBranch.SkeletonIDs)-1):-1:2
                            point_1 = tempBranch.Curve(m,:);
                            point_2 = tempBranch.Curve(m-1,:);
                            branchSegments_2 = [branchSegments_2;[point_1,point_2,associatedBranchIDs(l)]];
                        end
                    end
                end
                branchSegments_3 = [];
                for l = 1:1:length(associatedBranchIDs)
                    tempBranch = treeModel(associatedBranchIDs(l));
                    selectID = find(tempBranch.SkeletonIDs == bifurcationID);
                    if (associatedBranchIDs(l) ~= associatedBranchIDs(i)) && (selectID == 1)
                        point_1 = tempBranch.Curve(2,:);
                        point_2 = targetPoint;
                        branchSegments_3 = [branchSegments_3;[point_1,point_2,associatedBranchIDs(l)]];
                    end
                    if (associatedBranchIDs(l) ~= associatedBranchIDs(i)) && (selectID ~= 1)
                        point_1 = tempBranch.Curve(end-1,:);
                        point_2 = targetPoint;
                        branchSegments_3 = [branchSegments_3;[point_1,point_2,associatedBranchIDs(l)]];
                    end
                end
                branchSegments = [branchSegments_2;branchSegments_3];
                % Determine whether the adjustment credible
                isCredible = L1TreeIsBifurcationCredible(branchSegments,parentBranchID,targetPoint,parameters);
                if ~isCredible
                    continue;
                end
                newDistanceToTreeModel_1 = ones(size(associatedTLSPoints,1),1)*9999;
                totalLength_1 = 0;
                for l = 1:1:size(branchSegments_1,1)
                    point_1 = branchSegments_1(l,1:3);
                    point_2 = branchSegments_1(l,4:6);
                    P1P2 = point_1 - point_2;
                    P2P1 = point_2 - point_1;
                    P1P0 = point_1 - associatedTLSPoints;
                    tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                    tMatrix(tMatrix < 0) = 0;
                    tMatrix(tMatrix > 1) = 1;
                    pedalPoints = point_1 + tMatrix*P2P1;
                    pedalDistances = sqrt(sum((pedalPoints - associatedTLSPoints).^2,2));
                    selectIDs = find(pedalDistances < newDistanceToTreeModel_1);
                    newDistanceToTreeModel_1(selectIDs) = pedalDistances(selectIDs);
                    totalLength_1 = totalLength_1 + sqrt(sum((point_1 - point_2).^2));
                end
                newDistanceToTreeModel_2 = ones(size(associatedTLSPoints,1),1)*9999;
                totalLength_2 = 0;
                for l = 1:1:size(branchSegments_2,1)
                    point_1 = branchSegments_2(l,1:3);
                    point_2 = branchSegments_2(l,4:6);
                    P1P2 = point_1 - point_2;
                    P2P1 = point_2 - point_1;
                    P1P0 = point_1 - associatedTLSPoints;
                    tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                    tMatrix(tMatrix < 0) = 0;
                    tMatrix(tMatrix > 1) = 1;
                    pedalPoints = point_1 + tMatrix*P2P1;
                    pedalDistances = sqrt(sum((pedalPoints - associatedTLSPoints).^2,2));
                    selectIDs = find(pedalDistances < newDistanceToTreeModel_2);
                    newDistanceToTreeModel_2(selectIDs) = pedalDistances(selectIDs);
                    totalLength_2 = totalLength_2 + sqrt(sum((point_1 - point_2).^2));
                end
                newDistanceToTreeModel_3 = ones(size(associatedTLSPoints,1),1)*9999;
                totalLength_3 = 0;
                for l = 1:1:size(branchSegments_3,1)
                    point_1 = branchSegments_3(l,1:3);
                    point_2 = branchSegments_3(l,4:6);
                    P1P2 = point_1 - point_2;
                    P2P1 = point_2 - point_1;
                    P1P0 = point_1 - associatedTLSPoints;
                    tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                    tMatrix(tMatrix < 0) = 0;
                    tMatrix(tMatrix > 1) = 1;
                    pedalPoints = point_1 + tMatrix*P2P1;
                    pedalDistances = sqrt(sum((pedalPoints - associatedTLSPoints).^2,2));
                    selectIDs = find(pedalDistances < newDistanceToTreeModel_3);
                    newDistanceToTreeModel_3(selectIDs) = pedalDistances(selectIDs);
                    totalLength_3 = totalLength_3 + sqrt(sum((point_1 - point_2).^2));
                end
                associatedTLSPointIDs_1 = find(newDistanceToTreeModel_1 < newDistanceToTreeModel_2);
                associatedTLSPointIDs_2 = find(newDistanceToTreeModel_3 < newDistanceToTreeModel_2);
                associatedTLSPointIDs = unique([associatedTLSPointIDs_1;associatedTLSPointIDs_2]);
                meanDistances_1 = mean(newDistanceToTreeModel_1(associatedTLSPointIDs));
                meanDistances_3 = mean(newDistanceToTreeModel_3(associatedTLSPointIDs));
                if (abs((meanDistances_1 - meanDistances_3)/meanDistances_3) > abs((totalLength_3 - totalLength_1)/totalLength_1)) && (meanDistances_3 < meanDistances_1)
                    adjustmentStrategies = [adjustmentStrategies;[bifurcationID,bifurcationPoint,targetID,targetPoint,totalLength_2+totalLength_3,...
                    mean(min(newDistanceToTreeModel_1,newDistanceToTreeModel_3)),isTerminal,isInterpolation]];
                end
            end
        end
    end
    % Step 1.4: Identify the optimial adjustment strategy
    if isempty(adjustmentStrategies)
        bifurcationID = bifurcationID;
        bifurcationPoint = bifurcationPoint;
        targetID = bifurcationID;
        targetPoint = bifurcationPoint;
        isTerminal = false;
        isInterpolation = false;
    else
        selectID = find(adjustmentStrategies(:,10) == min(adjustmentStrategies(:,10)));
        bifurcationID = adjustmentStrategies(selectID,1);
        bifurcationPoint = adjustmentStrategies(selectID,2:4);
        targetID = adjustmentStrategies(selectID,5);
        targetPoint = adjustmentStrategies(selectID,6:8);
        isTerminal = adjustmentStrategies(selectID,11);
        isInterpolation = adjustmentStrategies(selectID,12);
    end
    if ~isTerminal
        newBifurcationPoints = [newBifurcationPoints;targetPoint];
    else
        delID = find((newBifurcationPoints(:,1) == targetPoint(1)) & (newBifurcationPoints(:,2) == targetPoint(2)) & (newBifurcationPoints(:,3) == targetPoint(3)));
        newBifurcationPoints(delID,:) = [];
    end
    % Step 1.5: Apply the optimal adjustment strategy
    if bifurcationID == targetID
        continue;
    end
    initialNeighborMatrix = [];
    for i = 1:1:size(treeModel,2)
        tempBranch = treeModel(i);
        if ~ismember(i,associatedBranchIDs)
            for j = 1:1:length(tempBranch.SkeletonIDs)-1
                point_1 = tempBranch.Curve(j,:);
                point_2 = tempBranch.Curve(j+1,:);
                distance = sqrt(sum((point_1 - point_2).^2));
                initialNeighborMatrix = [initialNeighborMatrix;[point_1,point_2,distance]];
            end
        else
            selectID_1 = find(tempBranch.SkeletonIDs == bifurcationID);
            selectID_2 = find((tempBranch.SkeletonIDs == targetID));
            if ~isempty(selectID_2)
                if isTerminal
                    continue;
                end
                if isInterpolation == 1
                    point_1 = targetPoint;
                    point_2 = tempBranch.Curve(selectID_2,:);
                    distance = sqrt(sum((point_1 - point_2).^2));
                    initialNeighborMatrix = [initialNeighborMatrix;[point_1,point_2,distance]];
                end
                if selectID_1 == 1
                    for j = selectID_2:1:length(tempBranch.SkeletonIDs)-1
                        point_1 = tempBranch.Curve(j,:);
                        point_2 = tempBranch.Curve(j+1,:);
                        distance = sqrt(sum((point_1 - point_2).^2));
                        initialNeighborMatrix = [initialNeighborMatrix;[point_1,point_2,distance]];
                    end
                else
                    for j = selectID_2:-1:2
                        point_1 = tempBranch.Curve(j,:);
                        point_2 = tempBranch.Curve(j-1,:);
                        distance = sqrt(sum((point_1 - point_2).^2));
                        initialNeighborMatrix = [initialNeighborMatrix;[point_1,point_2,distance]];
                    end
                end
            else
                if selectID_1 == 1
                    for j = length(tempBranch.SkeletonIDs):-1:2
                        point_1 = tempBranch.Curve(j,:);
                        if j == 2
                            point_2 = targetPoint;
                        else
                            point_2 = tempBranch.Curve(j-1,:);
                        end
                        distance = sqrt(sum((point_1 - point_2).^2));
                        initialNeighborMatrix = [initialNeighborMatrix;[point_1,point_2,distance]];
                    end
                else
                    for j = 1:1:length(tempBranch.SkeletonIDs)-1
                        point_1 = tempBranch.Curve(j,:);
                        if j == length(tempBranch.SkeletonIDs)-1
                            point_2 = targetPoint;
                        else
                            point_2 = tempBranch.Curve(j+1,:);
                        end
                        distance = sqrt(sum((point_1 - point_2).^2));
                        initialNeighborMatrix = [initialNeighborMatrix;[point_1,point_2,distance]];
                    end
                end
            end
        end
    end
    skeletonPoints = [];
    for i = 1:1:size(initialNeighborMatrix,1)
        point_1 = initialNeighborMatrix(i,1:3);
        point_2 = initialNeighborMatrix(i,4:6);
        if isempty(skeletonPoints)
            skeletonPoints = [skeletonPoints;point_1;point_2];
        else
            selectID = find((skeletonPoints(:,1) == point_1(1)) & (skeletonPoints(:,2) == point_1(2)) & (skeletonPoints(:,3) == point_1(3)));
            if isempty(selectID)
                skeletonPoints = [skeletonPoints;point_1];
            end
            selectID = find((skeletonPoints(:,1) == point_2(1)) & (skeletonPoints(:,2) == point_2(2)) & (skeletonPoints(:,3) == point_2(3)));
            if isempty(selectID)
                skeletonPoints = [skeletonPoints;point_2];
            end
        end
    end
    neighborMatrix = [];
    for i = 1:1:size(initialNeighborMatrix,1)
        point_1 = initialNeighborMatrix(i,1:3);
        point_2 = initialNeighborMatrix(i,4:6);
        pointID_1 = find((skeletonPoints(:,1) == point_1(1)) & (skeletonPoints(:,2) == point_1(2)) & (skeletonPoints(:,3) == point_1(3)));
        pointID_2 = find((skeletonPoints(:,1) == point_2(1)) & (skeletonPoints(:,2) == point_2(2)) & (skeletonPoints(:,3) == point_2(3)));
        distance = initialNeighborMatrix(i,7);
        neighborMatrix = [neighborMatrix;[pointID_1,pointID_2,distance]];
    end
    basePoint = treeModel(1).Curve(1,:);
    rootID = find((skeletonPoints(:,1) == basePoint(1)) & (skeletonPoints(:,2) == basePoint(2)) & (skeletonPoints(:,3) == basePoint(3)));
    newTreeModel = L1TreeUpdate(neighborMatrix,skeletonPoints,rootID);
    newTLSPoints = tlsPoints;
    for i = 1:1:size(treeModel,2)
        oriBranchID = i;
        oriStartPoint = treeModel(i).Curve(1,:);
        oriEndPoint = treeModel(i).Curve(end,:);
        isMatched = false;
        for j = 1:1:size(newTreeModel,2)
            newBranchID = j;
            newStartPoint = newTreeModel(j).Curve(1,:);
            newEndPoint = newTreeModel(j).Curve(end,:);
            if all(oriStartPoint == newStartPoint) && all(oriEndPoint == newEndPoint)
                isMatched = true;
                selectIDs = find(tlsPoints(:,4) == oriBranchID);
                newTLSPoints(selectIDs,4) = newBranchID;
                break;
            end
        end
        if ~isMatched
            selectIDs = find(tlsPoints(:,4) == oriBranchID);
            newTLSPoints(selectIDs,4) = 9999;
        end
    end
    associatedTLSPointIDs = find(newTLSPoints(:,4) == 9999);
    currentTLSPoints = tlsPoints(associatedTLSPointIDs,1:3);
    tlsPointsToBranches = ones(size(currentTLSPoints,1),2)*9999;
    for i = 1:1:size(newTreeModel,2)
        for j = 1:1:size(newTreeModel(i).Curve,1)-1
            point_1 = newTreeModel(i).Curve(j,:);
            point_2 = newTreeModel(i).Curve(j+1,:);
            P1P2 = point_1 - point_2;
            P2P1 = point_2 - point_1;
            P1P0 = point_1 - currentTLSPoints;
            tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
            tMatrix(tMatrix > 1) = 1;
            tMatrix(tMatrix < 0) = 0;
            pedalPoints = point_1 + tMatrix*P2P1;
            pedalDistances = sqrt(sum((pedalPoints - currentTLSPoints).^2,2));
            selectIDs = find(pedalDistances < tlsPointsToBranches(:,2));
            tlsPointsToBranches(selectIDs,1) = i;
            tlsPointsToBranches(selectIDs,2) = pedalDistances(selectIDs);
        end
    end
    newTLSPoints(associatedTLSPointIDs,4) = tlsPointsToBranches(:,1);
    treeModel = newTreeModel;
    tlsPoints = newTLSPoints;
    waitbar(2/5,globalBar,num2str(size(bifurcationInfo,1)));
end
selectIDs = find(tlsPoints(:,4) ~= 0);
[treeModel,removeBranches] = L1TreeRemoveBranches(treeModel,tlsPoints(selectIDs,1:3),parameters);

%% Step 2: Identify associated and unassociated TLS points
tlsPointsToBranches = ones(size(tlsPoints,1),2)*9999;
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
        selectIDs = find(pedalDistances < tlsPointsToBranches(:,2));
        tlsPointsToBranches(selectIDs,1) = i;
        tlsPointsToBranches(selectIDs,2) = pedalDistances(selectIDs);
    end
end
tlsPointsToRemovedBranches = ones(size(tlsPoints,1),1)*9999;
for i = 1:1:size(removeBranches,2)
    for j = 1:1:size(removeBranches(i).Curve,1)-1
        point_1 = removeBranches(i).Curve(j,:);
        point_2 = removeBranches(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - tlsPoints(:,1:3);
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - tlsPoints(:,1:3)).^2,2));
        selectIDs = find(pedalDistances < tlsPointsToRemovedBranches);
        tlsPointsToRemovedBranches(selectIDs) = pedalDistances(selectIDs);
    end
end
if isempty(removeBranches)
    selectIDs = find(tlsPoints(:,4) == 0);
    tlsPoints(:,4) = tlsPointsToBranches(:,1);
    tlsPoints(selectIDs,4) = 0;
else
    selectID_1 = find(tlsPointsToRemovedBranches < tlsPointsToBranches(:,2));
    selectID_2 = find(tlsPoints(:,4) == 0);
    tlsPoints(:,4) = tlsPointsToBranches(:,1);
    tlsPoints(selectID_1,4) = 0;
    tlsPoints(selectID_2,4) = 0;
end

end