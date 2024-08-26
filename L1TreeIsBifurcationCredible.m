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
% L1TreeIsBifurcationCredible.m    The function for determining whether a bifurcation adjustment strategy is credible or not
%
% Version 1.0
% Latest update     07 Feb 2024
% 
% INPUTS:
% branchSegments    Branch segments after applying the current adjustment strategy
% parentBranchID    The ID of parent branch
% targetPoint    The coordinates of the target point
% parameters    The parameter set
% 
% OUTPUTS:
% isCredible    Whether the adjustment strategy is credible or not
% ------------------------------------------------------------------------------

function isCredible = L1TreeIsBifurcationCredible(branchSegments,parentBranchID,targetPoint,parameters)

%% Identify key direction vectors
parentBranchSegments = branchSegments(branchSegments(:,7) == parentBranchID,:);
if isempty(parentBranchSegments)
    parentDirection = [];
else
    selectID_1 = find((parentBranchSegments(:,1) == targetPoint(1)) & (parentBranchSegments(:,2) == targetPoint(2)) & (parentBranchSegments(:,3) == targetPoint(3)));
    selectID_2 = find((parentBranchSegments(:,4) == targetPoint(1)) & (parentBranchSegments(:,5) == targetPoint(2)) & (parentBranchSegments(:,6) == targetPoint(3)));
    if ~isempty(selectID_1)
        parentDirection = parentBranchSegments(selectID_1,1:3) - parentBranchSegments(selectID_1,4:6);
    end
    if ~isempty(selectID_2)
        parentDirection = parentBranchSegments(selectID_2,4:6) - parentBranchSegments(selectID_2,1:3);
    end
    parentDirection = parentDirection/sqrt(sum(parentDirection.^2));
end
childBranchSegments = branchSegments(branchSegments(:,7) ~= parentBranchID,:);
childBranchIDs = unique(childBranchSegments(:,7));
childDirections_1 = [];
childDirections_2 = [];
for i = 1:1:length(childBranchIDs)
    selectID_1 = find((childBranchSegments(:,1) == targetPoint(1)) & (childBranchSegments(:,2) == targetPoint(2)) & (childBranchSegments(:,3) == targetPoint(3)) & (childBranchSegments(:,7) == childBranchIDs(i)));
    selectID_2 = find((childBranchSegments(:,4) == targetPoint(1)) & (childBranchSegments(:,5) == targetPoint(2)) & (childBranchSegments(:,6) == targetPoint(3)) & (childBranchSegments(:,7) == childBranchIDs(i)));
    if ~isempty(selectID_1)
        middlePoint = childBranchSegments(selectID_1,4:6);
        childDirection = childBranchSegments(selectID_1,4:6) - childBranchSegments(selectID_1,1:3);
        childDirection = childDirection/sqrt(sum(childDirection.^2));
        childDirections_1 = [childDirections_1;[childBranchIDs(i),childDirection]];
    end
    if ~isempty(selectID_2)
        middlePoint = childBranchSegments(selectID_2,1:3);
        childDirection = childBranchSegments(selectID_2,1:3) - childBranchSegments(selectID_2,4:6);
        childDirection = childDirection/sqrt(sum(childDirection.^2));
        childDirections_1 = [childDirections_1;[childBranchIDs(i),childDirection]];
    end
    selectID_1 = find((childBranchSegments(:,1) == middlePoint(1)) & (childBranchSegments(:,2) == middlePoint(2)) & (childBranchSegments(:,3) == middlePoint(3)) & (childBranchSegments(:,4) ~= targetPoint(1)));
    selectID_2 = find((childBranchSegments(:,4) == middlePoint(1)) & (childBranchSegments(:,5) == middlePoint(2)) & (childBranchSegments(:,6) == middlePoint(3)) & (childBranchSegments(:,1) ~= targetPoint(1)));
    if ~isempty(selectID_1)
        childDirection = childBranchSegments(selectID_1,4:6) - childBranchSegments(selectID_1,1:3);
        childDirection = childDirection/sqrt(sum(childDirection.^2));
        childDirections_2 = [childDirections_2;[childBranchIDs(i),childDirection]];
    end
    if ~isempty(selectID_2)
        childDirection = childBranchSegments(selectID_2,1:3) - childBranchSegments(selectID_2,4:6);
        childDirection = childDirection/sqrt(sum(childDirection.^2));
        childDirections_2 = [childDirections_2;[childBranchIDs(i),childDirection]];
    end
end

%% Calculate angles
angles_1 = [];
if ~isempty(parentDirection)
    for i = 1:1:size(childDirections_1,1)
        angle = acos(parentDirection*childDirections_1(i,2:4)')*180/pi;
        angles_1 = [angles_1,angle];
    end
end
for i = 1:1:size(childDirections_1,1)
    childBranchID = childDirections_1(i,1);
    if isempty(childDirections_2)
        continue;
    end
    selectID = find(childDirections_2(:,1) == childBranchID);
    if ~isempty(selectID)
        angle = acos(childDirections_1(i,2:4)*childDirections_2(selectID,2:4)')*180/pi;
        angles_1 = [angles_1,angle];
    end
end
angles_2 = [];
for i = 1:1:size(childDirections_1,1)-1
    for j = (i+1):1:size(childDirections_1,1)
        angle = acos(childDirections_1(i,2:4)*childDirections_1(j,2:4)')*180/pi;
        angles_2 = [angles_2,angle];
    end
end

%% Identify whether the current adjustment strategy is credible
if isempty(angles_1)
    isCredible_1 = true;
else
    isCredible_1 = all(angles_1 <= parameters.maxAngleThreshold);
end
if isempty(angles_2)
    isCredible_2 = true;
else
    isCredible_2 = all(angles_2 > parameters.minAngleThreshold);
end
isCredible = isCredible_1 & isCredible_2;

end