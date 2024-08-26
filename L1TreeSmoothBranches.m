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
% L1TreeSmoothBranches.m    The function for smoothing branches
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

function [treeModel,tlsPoints] = L1TreeSmoothBranches(treeModel,tlsPoints,parameters,globalBar)

waitbar(3/5,globalBar,'Smooth branches');

%% Step 1: Remove anomalous skelton points present in each branch
for i = 1:1:size(treeModel,2)
    newCurve = round(treeModel(i).Curve,6);
    if size(treeModel(i).Curve,6)
        treeModel(i).Curve = newCurve;
        continue;
    end
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
            angles = [angles,acos(direction_1*direction_2')*180/pi];
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

%% Step 2: Update the ID of the branch associated with each TLS point
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
                selectIDs_1 = find((tMatrix ~= 1) & (pedalDistances < tlsPointsToBranches(:,3)));
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
            selectIDs = find(tlsPoints(:,4) == j);
            searchRange = quantile(tlsPointsToBranches(selectIDs,3),0.9)*2;
        else
            selectIDs = find(tlsPoints(:,4) == j);
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
tlsPoints(:,4) = 9999;
for i = 1:1:size(treeModel,2)
    selectIDs = find((tlsPointsToBranches(:,1) == i) & (tlsPointsToBranches(:,2) == 0));
    tlsPoints(selectIDs,4) = i;
end
selectIDs = find(tlsPoints(:,4) == 9999);
tlsPoints(selectIDs,4) = 0;
associatedTLSPointIDs = sort(unique(associatedTLSPointIDs));
selectIDs = setdiff(1:size(tlsPoints,1),associatedTLSPointIDs);
tlsPoints(selectIDs,4) = 0;
[tlsPointIDs,distances] = rangesearch(tlsPoints(:,1:3),tlsPoints(:,1:3),parameters.initialSearchRange*2,'Distance','euclidean','NSMethod','kdtree');
tlsPoints = [tlsPoints,ones(size(tlsPoints,1),1)*9999];
for i = 1:1:size(tlsPoints,1)
    TLSNeighborIDs = setdiff(tlsPointIDs{i},i);
    if isempty(TLSNeighborIDs)
        tlsPoints(i,5) = 1;
    else
        tempDistances = setdiff(distances{i},0);
        weights = exp(tempDistances.^2*(-1/(parameters.initialSearchRange*2/2)^2));
        tlsPoints(i,5) = 1/(1+sum(weights));
    end
end

%% Step 3: Smooth branches
for i = 1:1:size(treeModel,2)
    selectIDs = find(tlsPoints(:,4) == i);
    currentTLSPoints = tlsPoints(selectIDs,1:3);
    currentTLSDensties = tlsPoints(selectIDs,5);
    tlsPointsToBranches = ones(size(currentTLSPoints,1),2)*9999;
    for j = 1:1:size(treeModel(i).Curve,1)-1
        point_1 = treeModel(i).Curve(j,:);
        point_2 = treeModel(i).Curve(j+1,:);
        P1P2 = point_1 - point_2;
        P2P1 = point_2 - point_1;
        P1P0 = point_1 - currentTLSPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix < 0) = 0;
        tMatrix(tMatrix > 1) = 1;
        pedalPoints = point_1 + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - currentTLSPoints).^2,2));
        selectIDs = find(pedalDistances < tlsPointsToBranches(:,2));
        if j == 1
            selectIDs_1 = find((tMatrix == 0) & (pedalDistances < tlsPointsToBranches(:,2)));
            selectIDs_2 = find((tMatrix ~= 0) & (pedalDistances < tlsPointsToBranches(:,2)));
            tlsPointsToBranches(selectIDs_1,1) = 1;
            tlsPointsToBranches(selectIDs_2,1) = 0;
        elseif j == size(treeModel(i).Curve,1)-1
            selectIDs_1 = find((tMatrix == 1) & (pedalDistances < tlsPointsToBranches(:,2)));
            selectIDs_2 = find((tMatrix ~= 1) & (pedalDistances < tlsPointsToBranches(:,2)));
            tlsPointsToBranches(selectIDs_1,1) = 1;
            tlsPointsToBranches(selectIDs_2,1) = 0;
        else
            tlsPointsToBranches(selectIDs,1) = 0;
        end
        tlsPointsToBranches(selectIDs,2) = pedalDistances(selectIDs);
    end
    selectIDs = find(tlsPointsToBranches(:,1) == 0);
    currentTLSPoints = currentTLSPoints(selectIDs,:);
    startPoint = treeModel(i).Curve(1,:);
    endPoint = treeModel(i).Curve(end,:);
    direction = endPoint - startPoint;
    direction = direction/sqrt(sum(direction.^2));
    distanceToStart = (currentTLSPoints - startPoint)*direction';
    [~,distanceOrder] = sort(distanceToStart);
    currentTLSPoints = currentTLSPoints(distanceOrder,:);
    currentTLSDensties = currentTLSDensties(distanceOrder,:);
    % Calculate the coordinates of new skeleton points
    if size(currentTLSPoints,1) < 10
        continue;
    end
    branchLength = sqrt(sum((endPoint - startPoint).^2));
    maxSampleNum = min(floor(branchLength/parameters.resolution),size(currentTLSPoints,1)/2);
    newCurve = [];
    tempRadius = 9999;
    for j = 0:1:maxSampleNum
        if j == 0
            tempCurve = [startPoint;endPoint];
        else
            tempCurve = startPoint;
            for k = 1:1:j
                startID = floor(size(currentTLSPoints,1)/(j+1)*(k-1))+1;
                endID = floor(size(currentTLSPoints,1)/(j+1)*(k+1));
                tempTLSPoints = currentTLSPoints(startID:endID,:);
                tempTLSDensities = currentTLSDensties(startID:endID,:);
                tempCurve = [tempCurve;round(tempTLSDensities'*tempTLSPoints/sum(tempTLSDensities),6)];
            end
            tempCurve = [tempCurve;endPoint];
        end
        isCredible = false;
        if size(tempCurve,1) < 3
            isCredible = true;
        else
            angles = [];
            for k = 2:1:size(tempCurve,1)-1
                direction_1 = tempCurve(k,:) - tempCurve(k-1,:);
                direction_1 = direction_1/sqrt(sum(direction_1.^2));
                direction_2 = tempCurve(k+1,:) - tempCurve(k,:);
                direction_2 = direction_2/sqrt(sum(direction_2.^2));
                angles = [angles,acos(direction_1*direction_2')*180/pi];
            end
            if max(angles) < parameters.smoothAngleThreshold
                isCredible = true;
            end
        end
        if isCredible
            tlsPointsToBranches = ones(size(currentTLSPoints,1),1)*9999;
            for k = 1:1:size(tempCurve,1)-1
                point_1 = tempCurve(k,:);
                point_2 = tempCurve(k+1,:);
                P1P2 = point_1 - point_2;
                P2P1 = point_2 - point_1;
                P1P0 = point_1 - currentTLSPoints;
                tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                tMatrix(tMatrix < 0) = 0;
                tMatrix(tMatrix > 1) = 1;
                pedalPoints = point_1 + tMatrix*P2P1;
                pedalDistances = sqrt(sum((pedalPoints - currentTLSPoints).^2,2));
                selectIDs = find(pedalDistances < tlsPointsToBranches);
                tlsPointsToBranches(selectIDs) = pedalDistances(selectIDs);
            end
        end
        if mean(tlsPointsToBranches) < tempRadius
            newCurve = tempCurve;
            tempRadius = mean(tlsPointsToBranches);
        end
    end
    if isempty(newCurve)
        newCurve = treeModel(i).Curve;
    end
    treeModel(i).Curve = newCurve;
end

%% Step 4: Update the 3D tree model
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

end