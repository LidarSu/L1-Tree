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
% L1MedianSearchNewBranches.m    The function for searching new branches
% 
% Version 1.0
% Latest update     28 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
% 
% INPUTS:
% skeletonPoints    The struct that records the attributes of skeleton points
% parameters    The parameter set
% globalBar    The global waiting bar
% iterationNum    The number of iterations
% maxIterationNum    The maximum number of iterations
% 
% OUTPUTS:
% newBranches    The struct that records the attributes of branches (the initial skeleton).
%                The attributes include the IDs of the skeleton points that make up each branch (SkeletonIDs);
%                and the coordinates of the skeleton points that make up each branch (Curve)
% skeletonPoints    The updated struct that records the attributes of skeleton points
% ------------------------------------------------------------------------------

function [newBranches,skeletonPoints] = L1MedianSearchNewBranches(skeletonPoints,parameters,globalBar,iterationNum,maxIterationNum)

waitbar(iterationNum/maxIterationNum,globalBar,'Search new branches');

%% Update skeletonPoints.NeighborSkeletonIDs and skeletonPoints.NeighborSkeletonPoints
skeletonPointPs = [skeletonPoints(:).P];
skeletonPointPs = reshape(skeletonPointPs,3,size(skeletonPoints,2))';
[~,kIDs] = pdist2(skeletonPointPs,skeletonPointPs,'euc','Smallest',parameters.kTracingThreshold+1);
kIDs = kIDs(2:end,:);
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsBranch || skeletonPoints(i).IsIgnore
        continue;
    end
    skeletonPoints(i).NeighborSkeletonIDs = [skeletonPoints(kIDs(:,i)).ID];
    tempNeighborSkeletonPoints = [skeletonPoints(kIDs(:,i)).P];
    skeletonPoints(i).NeighborSkeletonPoints = reshape(tempNeighborSkeletonPoints,3,length(skeletonPoints(i).NeighborSkeletonIDs))';
end

%% Search new branches
newBranches = [];
while true
    % Step 1: Identify the first point of a potential skeleton segment
    maxSigma = 0;
    firstID = -1;
    for i = 1:1:size(skeletonPoints,2)
        if skeletonPoints(i).IsFix && ~skeletonPoints(i).IsIgnore && ~skeletonPoints(i).IsBranch
            if skeletonPoints(i).Sigma > maxSigma
                maxSigma = skeletonPoints(i).Sigma;
                firstID = i;
            end
        end
    end
    if firstID == -1
        break;
    end
    tempNeighborSkeletonIDs = skeletonPoints(firstID).NeighborSkeletonIDs;
    secondID = -1;
    for i = 1:1:length(tempNeighborSkeletonIDs)
        if skeletonPoints(tempNeighborSkeletonIDs(i)).IsFix && ~skeletonPoints(tempNeighborSkeletonIDs(i)).IsIgnore && ~skeletonPoints(tempNeighborSkeletonIDs(i)).IsBranch
            secondID = tempNeighborSkeletonIDs(i);
            break;
        end
    end
    % Step 2: Extend the branch in different directions
    if secondID == -1
        skeletonPoints(firstID).IsFix = false;
        skeletonPoints(firstID).IsBranch = false;
        continue;
    end
    firstPoint = skeletonPoints(firstID).P;
    secondPoint = skeletonPoints(secondID).P;
    searchDirection = (secondPoint - firstPoint)/sqrt(sum((secondPoint - firstPoint).^2));
    [newBranch_1,skeletonPoints] = L1MedianExtendBranchViaDirection(firstID,searchDirection,skeletonPoints,parameters);
    searchDirection = -searchDirection;
    [newBranch_2,skeletonPoints] = L1MedianExtendBranchViaDirection(firstID,searchDirection,skeletonPoints,parameters);
    newBranch.SkeletonIDs = [];
    for i = length(newBranch_2.SkeletonIDs):-1:2
        newBranch.SkeletonIDs = [newBranch.SkeletonIDs,newBranch_2.SkeletonIDs(i)];
    end
    for i = 1:1:length(newBranch_1.SkeletonIDs)
        newBranch.SkeletonIDs = [newBranch.SkeletonIDs,newBranch_1.SkeletonIDs(i)];
    end
    if length(newBranch.SkeletonIDs) < parameters.branchSizeThreshold
        skeletonPoints(firstID).IsFix = false;
        skeletonPoints(firstID).IsBranch = false;
        continue;
    end
    %% Step 3: Extend the branch from ends
    [newBranch,skeletonPoints] = L1MedianExtendBranchFromTail(newBranch,skeletonPoints,parameters);
    tempBranch = newBranch;
    for i = 1:1:length(tempBranch.SkeletonIDs)
        tempBranch.SkeletonIDs(i) = newBranch.SkeletonIDs(length(tempBranch.SkeletonIDs)+1-i);
    end
    newBranch = tempBranch;
    [newBranch,skeletonPoints] = L1MedianExtendBranchFromTail(newBranch,skeletonPoints,parameters);
    tempBranch = newBranch;
    for i = 1:1:length(tempBranch.SkeletonIDs)
        tempBranch.SkeletonIDs(i) = newBranch.SkeletonIDs(length(tempBranch.SkeletonIDs)+1-i);
    end
    newBranch = tempBranch;
    %% Step 4: Update skeletonPoints.IsFix and skeletonPoints.IsBranch
    for i = 1:1:length(newBranch.SkeletonIDs)
        skeletonPoints(newBranch.SkeletonIDs(i)).IsFix = true;
        skeletonPoints(newBranch.SkeletonIDs(i)).IsBranch = true;
    end
    newBranch.Curve = zeros(length(newBranch.SkeletonIDs),3);
    for i = 1:1:length(newBranch.SkeletonIDs)
        newBranch.Curve(i,:) = skeletonPoints(newBranch.SkeletonIDs(i)).P;
    end
    newBranches = [newBranches,newBranch];
end

end