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
% L1TreeIsConnectionCredible.m    The function for determining whether a potential connection is credible or not
%
% Version 1.0
% Latest update     03 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% connection    The potential connection
% targetComponent    The target connected component
% allSkeletonPoints    The coordinates of all skeleton points
% tlsPoints    The coordinates of TLS points
% ignoreDirections    Vector directions that should be ignored when applying angle constraints
% tlsPointsToComponents    The table recording the position relationship between TLS points and connected components.
%                          The first column records the ID of the nearest connected component to a TLS point;
%                          the second column records the ID of the nearest skeleton point to a TLS point;
%                          the third column records whether the nearest skeleton point to a TLS point is an end point or not;
%                          and the forth column records the distance between a TLS point to its nearest connected component
% connectionParameters    The struct recording the information about the potential connection.
%                         The information includes the ID of the connected component under examination (CurrentComponentID);
%                         the connected component (CurrentComponent);
%                         the direction at the end point of the current connected component (startDirection);
%                         the IDs of TLS points associated with the last connected component (LastAddTLSPointIDs); 
%                         and the search range determined based on the last connected component (LastSearchRange)
% parameters    The parameter set
% 
% OUTPUTS:
% isCredible    Whether the connection is credible or not
% connectionStrategy    The struct recording the information about the connection strategy.
%                       The information includes the ID of the connected component under examination (CurrentComponentID);
%                       the potential connection (connection); the direction at the end point of the current connected component (startDirection);
%                       number of TLS points with a shorter distance to the connection than other branches (RelatedNum);
%                       the length of the connection (Length); the average distance between the connection and its associated TLS points;
%                       the IDs of TLS points associated with the last connected component (LastAddTLSPointIDs); 
%                       and the search range determined based on the last connected component (LastSearchRange)
% ------------------------------------------------------------------------------

function [isCredible,connectionStrategy] = L1TreeIsConnectionCredible(connection,targetComponent,allSkeletonPoints,tlsPoints,ignoreDirections,tlsPointsToComponents,connectionParameters,parameters)

if length(connection) == 2
    % Step 1: Obtain basic information about the connection
    currentComponentEndPointID = connection(1);
    lastComponentPointID = connection(2);
    currentComponentEndPoint = allSkeletonPoints(currentComponentEndPointID,:);
    lastComponentPoint = allSkeletonPoints(lastComponentPointID,:);
    connectionDirection = lastComponentPoint - currentComponentEndPoint;
    connectionDirection = connectionDirection/sqrt(sum(connectionDirection.^2));
    currentComponentID = connectionParameters.CurrentComponentID;
    currentComponent = connectionParameters.CurrentComponent;
    startDirection = connectionParameters.StartDirection;
    lastSearchRange = connectionParameters.LastSearchRange;
    % Step 2: Identify three sets of directions which are crucial for evaluating the credible of the connection
    directions_1 = [];
    for i = 1:1:size(targetComponent,2)
        selectID = find(targetComponent(i).SkeletonIDs == lastComponentPointID);
        if isempty(selectID)
            continue;
        end
        if selectID == 1
            for j = 1:1:size(targetComponent,2)
                if contains(targetComponent(i).CenCode,targetComponent(j).CenCode) && (targetComponent(j).CenOrder == targetComponent(i).CenOrder-1)
                    pointID_1 = targetComponent(j).SkeletonIDs(1);
                    pointID_2 = targetComponent(j).SkeletonIDs(end-1);
                    pointID_3 = targetComponent(j).SkeletonIDs(end);
                    point_1 = allSkeletonPoints(pointID_1,:);
                    point_2 = allSkeletonPoints(pointID_2,:);
                    point_3 = allSkeletonPoints(pointID_3,:);
                    direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                    directions_1 = [directions_1;[pointID_1,pointID_3,direction]];
                    direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                    directions_1 = [directions_1;[pointID_2,pointID_3,direction]];
                    break;
                end
            end
        end
        if selectID == length(targetComponent(i).SkeletonIDs)
            pointID_1 = targetComponent(i).SkeletonIDs(1);
            pointID_2 = targetComponent(i).SkeletonIDs(end-1);
            pointID_3 = targetComponent(i).SkeletonIDs(end);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_1 = [directions_1;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_1 = [directions_1;[pointID_2,pointID_3,direction]];
        end
        if (selectID ~= 1) && (selectID ~= length(targetComponent(i).SkeletonIDs))
            pointID_1 = targetComponent(i).SkeletonIDs(1);
            pointID_2 = targetComponent(i).SkeletonIDs(selectID-1);
            pointID_3 = targetComponent(i).SkeletonIDs(selectID);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_1 = [directions_1;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_1 = [directions_1;[pointID_2,pointID_3,direction]];
        end
        break;
    end
    if ~isempty(directions_1)
        delIDs = find(ismember(directions_1(:,1),ignoreDirections(:,1)) & ismember(directions_1(:,2),ignoreDirections(:,2)));
        directions_1(delIDs,:) = [];
        directions_1 = directions_1(:,3:5);
    end
    directions_2 = [];
    for i = 1:1:size(targetComponent,2)
        selectID = find(targetComponent(i).SkeletonIDs == lastComponentPointID);
        if isempty(selectID)
            continue;
        end
        if selectID == 1
            parentBranchID = [];
            for j = 1:1:size(targetComponent,2)
                if contains(targetComponent(i).CenCode,targetComponent(j).CenCode) && (targetComponent(j).CenCode == targetComponent(j).CenOrder-1)
                    parentBranchID = j;
                    break;
                end
            end
            if isempty(parentBranchID)
                pointID_1 = targetComponent(i).SkeletonIDs(end);
                pointID_2 = targetComponent(i).SkeletonIDs(2);
                pointID_3 = targetComponent(i).SkeletonIDs(1);
                point_1 = allSkeletonPoints(pointID_1,:);
                point_2 = allSkeletonPoints(pointID_2,:);
                point_3 = allSkeletonPoints(pointID_3,:);
                direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
                direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
            else
                for j = 1:1:size(targetComponent,2)
                    if contains(targetComponent(j).CenCode,targetComponent(parentBranchID).CenCode) && (targetComponent(j).CenOrder == targetComponent(parentBranchID).CenOrder+1)
                        pointID_1 = targetComponent(j).SkeletonIDs(end);
                        pointID_2 = targetComponent(j).SkeletonIDs(2);
                        pointID_3 = targetComponent(j).SkeletonIDs(1);
                        point_1 = allSkeletonPoints(pointID_1,:);
                        point_2 = allSkeletonPoints(pointID_2,:);
                        point_3 = allSkeletonPoints(pointID_3,:);
                        direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                        directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
                        direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                        directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
                    end
                end
            end
        end
        if selectID == length(targetComponent(i).SkeletonIDs)
            for j = 1:1:size(targetComponent,2)
                if contains(targetComponent(j).CenCode,targetComponent(i).CenCode) && (targetComponent(j).CenOrder == targetComponent(i).CenOrder+1)
                    pointID_1 = targetComponent(j).SkeletonIDs(end);
                    pointID_2 = targetComponent(j).SkeletonIDs(2);
                    pointID_3 = targetComponent(j).SkeletonIDs(1);
                    point_1 = allSkeletonPoints(pointID_1,:);
                    point_2 = allSkeletonPoints(pointID_2,:);
                    point_3 = allSkeletonPoints(pointID_3,:);
                    direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                    directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
                    direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                    directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
                end
            end
        end
        if (selectID ~= 1) && (selectID ~= length(targetComponent(i).SkeletonIDs))
            pointID_1 = targetComponent(i).SkeletonIDs(end);
            pointID_2 = targetComponent(i).SkeletonIDs(selectID+1);
            pointID_3 = targetComponent(i).SkeletonIDs(selectID);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
        end
        break;
    end
    if ~isempty(directions_2)
        delIDs = find(ismember(directions_2(:,1),ignoreDirections(:,1)) & ismember(directions_2(:,2),ignoreDirections(:,2)));
        directions_2(delIDs,:) = [];
        directions_2 = directions_2(:,3:5);
    end
    directions_3 = [];
    for i = 1:1:size(currentComponent,2)
        selectID = find(currentComponent(i).SkeletonIDs == currentComponentEndPointID);
        if isempty(selectID)
            continue;
        end
        if selectID == 1
            pointID_1 = currentComponent(i).SkeletonIDs(end);
            pointID_2 = currentComponent(i).SkeletonIDs(2);
            pointID_3 = currentComponent(i).SkeletonIDs(1);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_3 = [directions_3;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_3 = [directions_3;[pointID_2,pointID_3,direction]];
        end
        if selectID == length(currentComponent(i).SkeletonIDs)
            pointID_1 = currentComponent(i).SkeletonIDs(1);
            pointID_2 = currentComponent(i).SkeletonIDs(end-1);
            pointID_3 = currentComponent(i).SkeletonIDs(end);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_3 = [directions_3;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_3 = [directions_3;[pointID_2,pointID_3,direction]];
        end
        if (selectID ~= 1) && (selectID ~= length(currentComponent(i).SkeletonIDs))
            pointID_1 = currentComponent(i).SkeletonIDs(1);
            pointID_2 = currentComponent(i).SkeletonIDs(selectID-1);
            pointID_3 = currentComponent(i).SkeletonIDs(selectID);
            pointID_4 = currentComponent(i).SkeletonIDs(selectID+1);
            pointID_5 = currentComponent(i).SkeletonIDs(end);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            point_4 = allSkeletonPoints(pointID_4,:);
            point_5 = allSkeletonPoints(pointID_5,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_3 = [directions_3;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_3 = [directions_3;[pointID_2,pointID_3,direction]];
            direction = (point_4 - point_3)/sqrt(sum((point_4 - point_3).^2));
            directions_3 = [directions_3;[pointID_4,pointID_3,direction]];
            direction = (point_5 - point_3)/sqrt(sum((point_5 - point_3).^2));
            directions_3 = [directions_3;[pointID_5,pointID_3,direction]];
        end
    end
    if ~isempty(directions_3)
        delIDs = find(ismember(directions_3(:,1),ignoreDirections(:,1)) & ismember(directions_3(:,2),ignoreDirections(:,2)));
        directions_3(delIDs,:) = [];
        directions_3 = directions_3(:,3:5);
    end
    % Step 3: Apply angle constraints
    angles_1 = [];
    angle = acos(startDirection*connectionDirection')*180/pi;
    angles_1 = [angles_1,angle];
    for i = 1:1:size(directions_1,1)
        angle = acos(directions_1(i,:)*connectionDirection')*180/pi;
        angles_1 = [angles_1,angle];
    end
    isCredible_1 = all(angles_1 <= parameters.maxAngleThreshold);
    angles_2 = [];
    for i = 1:1:size(directions_2,1)
        angle = acos(directions_2(i,:)*(-connectionDirection)')*180/pi;
        angles_2 = [angles_2,angle];
    end
    if isempty(angles_2)
        isCredible_2 = true;
    else
        isCredible_2 = all(angles_2 > parameters.minAngleThreshold);
    end
    angles_3 = [];
    for i = 1:1:size(directions_3,1)
        angle = acos(directions_3(i,:)*(-connectionDirection)')*180/pi;
        angles_3 = [angles_3,angle];
    end
    if isempty(angles_3)
        isCredible_3 = true;
    else
        isCredible_3 = all(angles_3 <= parameters.maxAngleThreshold);
    end
    angles_4 = [];
    for i = 1:1:size(directions_3,1)
        angle = acos(directions_3(i,:)*connectionDirection')*180/pi;
        angles_4 = [angles_4,angle];
    end
    if isempty(angles_4)
        isCredible_4 = true;
    else
        isCredible_4 = all(angles_4 > parameters.minAngleThreshold);
    end
    isCredible = isCredible_1 & isCredible_2 & isCredible_3 & isCredible_4;
    if ~isCredible
        connectionStrategy = [];
        return;
    end
    % Step 4: Apply point density constraints
    connectionLength = sqrt(sum((lastComponentPoint - currentComponentEndPoint).^2));
    P1P2 = lastComponentPoint - currentComponentEndPoint;
    P2P1 = currentComponentEndPoint - lastComponentPoint;
    P1P0 = lastComponentPoint - tlsPoints;
    tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
    tMatrix(tMatrix > 1) = 1;
    tMatrix(tMatrix < 0) = 0;
    pedalPoints = lastComponentPoint + tMatrix*P2P1;
    pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
    relatedTLSPointIDs_1 = find(tlsPointsToComponents(:,3) == 1);
    relatedTLSPointIDs_2 = connectionParameters.LastAddTLSPointIDs;
    relatedTLSPointIDs = find((pedalDistances < tlsPointsToComponents(:,4)) & (pedalDistances < lastSearchRange*2) & (tMatrix ~= 0) & (tMatrix ~= 1));
    relatedTLSPointIDs = relatedTLSPointIDs(ismember(relatedTLSPointIDs,[relatedTLSPointIDs_1;relatedTLSPointIDs_2]));
    pedalPoints = pedalPoints(relatedTLSPointIDs,:);
    pedalDistances = pedalDistances(relatedTLSPointIDs);
    if isempty(relatedTLSPointIDs)
        isCredible = false;
    else
        isCredible = true;
        if connectionLength > parameters.minDistanceThreshold
            sampleNum = floor((connectionLength - parameters.windowSize)/parameters.moveStep);
            for i = 0:1:sampleNum
                P1 = lastComponentPoint + (currentComponentEndPoint - lastComponentPoint)/connectionLength*i*parameters.moveStep;
                P2 = P1 + (currentComponentEndPoint - lastComponentPoint)/connectionLength*parameters.windowSize;
                if P1(1) > P2(1)
                    selectIDs = find((pedalPoints(:,1) > P2(1)) & (pedalPoints(:,1) < P1(1)));
                else
                    selectIDs = find((pedalPoints(:,1) > P1(1)) & (pedalPoints(:,1) < P2(1)));
                end
                if length(selectIDs) <= length(relatedTLSPointIDs)/connectionLength*parameters.windowSize*parameters.densityFactor
                    isCredible = false;
                    break;
                end
            end
        end
    end
    % Step 5: Output
    if ~isCredible
        connectionStrategy = [];
        return;
    else
        connectionStrategy.CurrentComponentID = currentComponentID;
        connectionStrategy.Connection = connection;
        connectionStrategy.StartDirection = startDirection;
        connectionStrategy.RelatedNum = length(relatedTLSPointIDs);
        connectionStrategy.Length = connectionLength;
        connectionStrategy.Distance = mean(pedalDistances);
        connectionStrategy.LastAddTLSPointIDs = relatedTLSPointIDs_2;
        connectionStrategy.LastSearchRange = lastSearchRange;
        return;
    end
end

if length(connection) > 2
    % Step 1: Obtain basic information about the connection
    currentComponentEndPointID = connection(1);
    lastComponentPointID = connection(end);
    connectionDirection_1 = allSkeletonPoints(connection(end),:) - allSkeletonPoints(connection(1),:);
    connectionDirection_1 = connectionDirection_1/sqrt(sum(connectionDirection_1.^2));
    connectionDirection_2 = allSkeletonPoints(connection(2),:) - allSkeletonPoints(connection(1),:);
    connectionDirection_2 = connectionDirection_2/sqrt(sum(connectionDirection_2.^2));
    connectionDirection_3 = allSkeletonPoints(connection(end),:) - allSkeletonPoints(connection(end-1),:);
    connectionDirection_3 = connectionDirection_3/sqrt(sum(connectionDirection_3.^2));
    connectionDirections_1 = [connectionDirection_1;connectionDirection_2];
    connectionDirections_2 = [connectionDirection_1;connectionDirection_3];
    currentComponentID = connectionParameters.CurrentComponentID;
    currentComponent = connectionParameters.CurrentComponent;
    startDirection = connectionParameters.StartDirection;
    lastSearchRange = connectionParameters.LastSearchRange;
    % Step 2: Identify three sets of directions which are crucial for evaluating the credible of the connection
    directions_1 = [];
    for i = 1:1:size(targetComponent,2)
        selectID = find(targetComponent(i).SkeletonIDs == lastComponentPointID);
        if isempty(selectID)
            continue;
        end
        if selectID == 1
            for j = 1:1:size(targetComponent,2)
                if contains(targetComponent(i).CenCode,targetComponent(j).CenCode) && (targetComponent(j).CenOrder == targetComponent(i).CenOrder-1)
                    pointID_1 = targetComponent(j).SkeletonIDs(1);
                    pointID_2 = targetComponent(j).SkeletonIDs(end-1);
                    pointID_3 = targetComponent(j).SkeletonIDs(end);
                    point_1 = allSkeletonPoints(pointID_1,:);
                    point_2 = allSkeletonPoints(pointID_2,:);
                    point_3 = allSkeletonPoints(pointID_3,:);
                    direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                    directions_1 = [directions_1;[pointID_1,pointID_3,direction]];
                    direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                    directions_1 = [directions_1;[pointID_2,pointID_3,direction]];
                    break;
                end
            end
        end
        if selectID == length(targetComponent(i).SkeletonIDs)
            pointID_1 = targetComponent(i).SkeletonIDs(1);
            pointID_2 = targetComponent(i).SkeletonIDs(end-1);
            pointID_3 = targetComponent(i).SkeletonIDs(end);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_1 = [directions_1;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_1 = [directions_1;[pointID_2,pointID_3,direction]];
        end
        if (selectID ~= 1) && (selectID ~= length(targetComponent(i).SkeletonIDs))
            pointID_1 = targetComponent(i).SkeletonIDs(1);
            pointID_2 = targetComponent(i).SkeletonIDs(selectID-1);
            pointID_3 = targetComponent(i).SkeletonIDs(selectID);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_1 = [directions_1;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_1 = [directions_1;[pointID_2,pointID_3,direction]];
        end
        break;
    end
    if ~isempty(directions_1)
        delIDs = find(ismember(directions_1(:,1),ignoreDirections(:,1)) & ismember(directions_1(:,2),ignoreDirections(:,2)));
        directions_1(delIDs,:) = [];
        directions_1 = directions_1(:,3:5);
    end
    directions_2 = [];
    for i = 1:1:size(targetComponent,2)
        selectID = find(targetComponent(i).SkeletonIDs == lastComponentPointID);
        if isempty(selectID)
            continue;
        end
        if selectID == 1
            parentBranchID = [];
            for j = 1:1:size(targetComponent,2)
                if contains(targetComponent(i).CenCode,targetComponent(j).CenCode) && (targetComponent(j).CenCode == targetComponent(j).CenOrder-1)
                    parentBranchID = j;
                    break;
                end
            end
            if isempty(parentBranchID)
                pointID_1 = targetComponent(i).SkeletonIDs(end);
                pointID_2 = targetComponent(i).SkeletonIDs(2);
                pointID_3 = targetComponent(i).SkeletonIDs(1);
                point_1 = allSkeletonPoints(pointID_1,:);
                point_2 = allSkeletonPoints(pointID_2,:);
                point_3 = allSkeletonPoints(pointID_3,:);
                direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
                direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
            else
                for j = 1:1:size(targetComponent,2)
                    if contains(targetComponent(j).CenCode,targetComponent(parentBranchID).CenCode) && (targetComponent(j).CenOrder == targetComponent(parentBranchID).CenOrder+1)
                        pointID_1 = targetComponent(j).SkeletonIDs(end);
                        pointID_2 = targetComponent(j).SkeletonIDs(2);
                        pointID_3 = targetComponent(j).SkeletonIDs(1);
                        point_1 = allSkeletonPoints(pointID_1,:);
                        point_2 = allSkeletonPoints(pointID_2,:);
                        point_3 = allSkeletonPoints(pointID_3,:);
                        direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                        directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
                        direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                        directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
                    end
                end
            end
        end
        if selectID == length(targetComponent(i).SkeletonIDs)
            for j = 1:1:size(targetComponent,2)
                if contains(targetComponent(j).CenCode,targetComponent(i).CenCode) && (targetComponent(j).CenOrder == targetComponent(i).CenOrder+1)
                    pointID_1 = targetComponent(j).SkeletonIDs(end);
                    pointID_2 = targetComponent(j).SkeletonIDs(2);
                    pointID_3 = targetComponent(j).SkeletonIDs(1);
                    point_1 = allSkeletonPoints(pointID_1,:);
                    point_2 = allSkeletonPoints(pointID_2,:);
                    point_3 = allSkeletonPoints(pointID_3,:);
                    direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
                    directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
                    direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
                    directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
                end
            end
        end
        if (selectID ~= 1) && (selectID ~= length(targetComponent(i).SkeletonIDs))
            pointID_1 = targetComponent(i).SkeletonIDs(end);
            pointID_2 = targetComponent(i).SkeletonIDs(selectID+1);
            pointID_3 = targetComponent(i).SkeletonIDs(selectID);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_2 = [directions_2;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_2 = [directions_2;[pointID_2,pointID_3,direction]];
        end
        break;
    end
    if ~isempty(directions_2)
        delIDs = find(ismember(directions_2(:,1),ignoreDirections(:,1)) & ismember(directions_2(:,2),ignoreDirections(:,2)));
        directions_2(delIDs,:) = [];
        directions_2 = directions_2(:,3:5);
    end
    directions_3 = [];
    for i = 1:1:size(currentComponent,2)
        selectID = find(currentComponent(i).SkeletonIDs == currentComponentEndPointID);
        if isempty(selectID)
            continue;
        end
        if selectID == 1
            pointID_1 = currentComponent(i).SkeletonIDs(end);
            pointID_2 = currentComponent(i).SkeletonIDs(2);
            pointID_3 = currentComponent(i).SkeletonIDs(1);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_3 = [directions_3;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_3 = [directions_3;[pointID_2,pointID_3,direction]];
        end
        if selectID == length(currentComponent(i).SkeletonIDs)
            pointID_1 = currentComponent(i).SkeletonIDs(1);
            pointID_2 = currentComponent(i).SkeletonIDs(end-1);
            pointID_3 = currentComponent(i).SkeletonIDs(end);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_3 = [directions_3;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_3 = [directions_3;[pointID_2,pointID_3,direction]];
        end
        if (selectID ~= 1) && (selectID ~= length(currentComponent(i).SkeletonIDs))
            pointID_1 = currentComponent(i).SkeletonIDs(1);
            pointID_2 = currentComponent(i).SkeletonIDs(selectID-1);
            pointID_3 = currentComponent(i).SkeletonIDs(selectID);
            pointID_4 = currentComponent(i).SkeletonIDs(selectID+1);
            pointID_5 = currentComponent(i).SkeletonIDs(end);
            point_1 = allSkeletonPoints(pointID_1,:);
            point_2 = allSkeletonPoints(pointID_2,:);
            point_3 = allSkeletonPoints(pointID_3,:);
            point_4 = allSkeletonPoints(pointID_4,:);
            point_5 = allSkeletonPoints(pointID_5,:);
            direction = (point_1 - point_3)/sqrt(sum((point_1 - point_3).^2));
            directions_3 = [directions_3;[pointID_1,pointID_3,direction]];
            direction = (point_2 - point_3)/sqrt(sum((point_2 - point_3).^2));
            directions_3 = [directions_3;[pointID_2,pointID_3,direction]];
            direction = (point_4 - point_3)/sqrt(sum((point_4 - point_3).^2));
            directions_3 = [directions_3;[pointID_4,pointID_3,direction]];
            direction = (point_5 - point_3)/sqrt(sum((point_5 - point_3).^2));
            directions_3 = [directions_3;[pointID_5,pointID_3,direction]];
        end
    end
    if ~isempty(directions_3)
        delIDs = find(ismember(directions_3(:,1),ignoreDirections(:,1)) & ismember(directions_3(:,2),ignoreDirections(:,2)));
        directions_3(delIDs,:) = [];
        directions_3 = directions_3(:,3:5);
    end
    % Step 3: Apply angle constraints
    angles_1 = [];
    for i = 1:1:size(connectionDirections_1,1)
        angle = acos(startDirection*connectionDirections_1(i,:)')*180/pi;
        angles_1 = [angles_1,angle];
    end
    for i = 1:1:size(directions_1,1)
        for j = 1:1:size(connectionDirections_2,1)
            angle = acos(directions_1(i,:)*connectionDirections_2(j,:)')*180/pi;
            angles_1 = [angles_1,angle];
        end
    end
    isCredible_1 = all(angles_1 <= parameters.maxAngleThreshold);
    angles_2 = [];
    for i = 1:1:size(directions_2,1)
        for j = 1:1:size(connectionDirections_2,1)
            angle = acos(directions_2(i,:)*(-connectionDirections_2(j,:))')*180/pi;
            angles_2 = [angles_2,angle];
        end
    end
    if isempty(angles_2)
        isCredible_2 = true;
    else
        isCredible_2 = all(angles_2 > parameters.minAngleThreshold);
    end
    angles_3 = [];
    for i = 1:1:size(directions_3,1)
        for j = 1:1:size(connectionDirections_1,1)
            angle = acos(directions_3(i,:)*(-connectionDirections_1(j,:))')*180/pi;
            angles_3 = [angles_3,angle];
        end
    end
    if isempty(angles_3)
        isCredible_3 = true;
    else
        isCredible_3 = all(angles_3 <= parameters.maxAngleThreshold);
    end
    angles_4 = [];
    for i = 1:1:size(directions_3,1)
        for j = 1:1:size(connectionDirections_1,1)
            angle = acos(directions_3(i,:)*connectionDirections_1(j,:)')*180/pi;
            angles_4 = [angles_4,angle];
        end
    end
    if isempty(angles_4)
        isCredible_4 = true;
    else
        isCredible_4 = all(angles_4 > parameters.minAngleThreshold);
    end
    isCredible = isCredible_1 & isCredible_2 & isCredible_3 & isCredible_4;
    if ~isCredible
        connectionStrategy = [];
        return;
    end
    % Step 4: Apply point density constraints
    totalLength = 0;
    relatedTLSPointIDs_1 = find(tlsPointsToComponents(:,3) == 1);
    relatedTLSPointIDs_2 = connectionParameters.LastAddTLSPointIDs;
    relatedTLSPointIDs = [];
    for i = 1:1:length(connection)-1
        P1P2 = allSkeletonPoints(connection(i+1),:) - allSkeletonPoints(connection(i),:);
        P2P1 = allSkeletonPoints(connection(i),:) - allSkeletonPoints(connection(i+1),:);
        P1P0 = allSkeletonPoints(connection(i+1),:) - tlsPoints;
        tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
        tMatrix(tMatrix > 1) = 1;
        tMatrix(tMatrix < 0) = 0;
        pedalPoints = allSkeletonPoints(connection(i+1),:) + tMatrix*P2P1;
        pedalDistances = sqrt(sum((pedalPoints - tlsPoints).^2,2));
        relatedTLSPointIDs = [relatedTLSPointIDs;find((pedalDistances < tlsPointsToComponents(:,4)) & (pedalDistances < lastSearchRange*2) & (tMatrix ~= 0) & (tMatrix ~= 1))];
        totalLength = totalLength + sqrt(sum(P1P2.^2));
    end
    relatedTLSPointIDs = unique(relatedTLSPointIDs);
    relatedTLSPointIDs = relatedTLSPointIDs(ismember(relatedTLSPointIDs,[relatedTLSPointIDs_1;relatedTLSPointIDs_2]));
    if isempty(relatedTLSPointIDs)
        isCredible = false;
    else
        isCredible = true;
        summaryTable = [];
        for i = 1:1:length(relatedTLSPointIDs)
            distance = 9999;
            for j = 1:1:length(connection)-1
                P1P2 = allSkeletonPoints(connection(j+1),:) - allSkeletonPoints(connection(j),:);
                P2P1 = allSkeletonPoints(connection(j),:) - allSkeletonPoints(connection(j+1),:);
                P1P0 = allSkeletonPoints(connection(j+1),:) - tlsPoints(relatedTLSPointIDs(i),:);
                tMatrix = (P1P0*P1P2')/(P1P2*P1P2');
                tMatrix(tMatrix > 1) = 1;
                tMatrix(tMatrix < 0) = 0;
                pedalPoint = allSkeletonPoints(connection(j+1),:) + tMatrix*P2P1;
                pedalDistance = sqrt(sum((pedalPoint - tlsPoints(relatedTLSPointIDs(i),:)).^2));
                if pedalDistance < distance
                    distance = pedalDistance;
                    tempTable = [j,pedalPoint,pedalDistance];
                end
            end
            if distance ~= 9999
                summaryTable = [summaryTable;tempTable];
            end
        end
        distances = [];
        for i = 1:1:length(relatedTLSPointIDs)
            beforeID = summaryTable(i,1);
            pedalPoint = summaryTable(i,2:4);
            distance = 0;
            if beforeID > 1
                for j = beforeID:-1:2
                    distance = distance + sqrt(sum((allSkeletonPoints(connection(j),:) - allSkeletonPoints(connection(j-1),:)).^2));
                end
            end
            distance = distance + sqrt(sum((pedalPoint - allSkeletonPoints(connection(beforeID),:)).^2));
            distances = [distances;distance];
        end
        summaryTable = [summaryTable,distances];
        if totalLength > parameters.minDistanceThreshold
            sampleNum = (totalLength - parameters.windowSize)/parameters.moveStep;
            for i = 0:1:sampleNum
                minDistance = i*parameters.moveStep;
                maxDistance = minDistance + parameters.windowSize;
                selectIDs = find((summaryTable(:,end) >= minDistance) & (summaryTable(:,end) < maxDistance));
                if length(selectIDs) <= length(relatedTLSPointIDs)/totalLength*parameters.windowSize*parameters.densityFactor
                    isCredible = false;
                    break;
                end
            end
        end
    end
    % Step 5: Output
    if ~isCredible
        connectionStrategy = [];
        return;
    else
        connectionStrategy.CurrentComponentID = currentComponentID;
        connectionStrategy.Connection = connection;
        connectionStrategy.StartDirection = startDirection;
        connectionStrategy.RelatedNum = length(relatedTLSPointIDs);
        connectionStrategy.Length = totalLength;
        connectionStrategy.Distance = mean(summaryTable(:,end-1));
        connectionStrategy.LastAddTLSPointIDs = relatedTLSPointIDs_2;
        connectionStrategy.LastSearchRange = lastSearchRange;
        return;
    end
end

end