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
% L1MedianExtendBranchViaDirection.m    Extend the branch in a specific direction
% 
% Version 1.0
% Latest update     28 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
% 
% INPUTS:
% firstID    The ID of the first skeleton point of a potential skeleton segment
% searchDirection    The search direction
% skeletonPoints    The struct that records the attributes of skeleton points
% parameters    The updated parameter set
% 
% OUTPUTS:
% newBranch    The obtained branch which includes the IDs of its constituent skeleton points
% skeletonPoints    The updated struct that records the attributes of skeleton points
% ------------------------------------------------------------------------------

function [newBranch,skeletonPoints] = L1MedianExtendBranchViaDirection(firstID,searchDirection,skeletonPoints,parameters)

currentID = firstID;
newBranch.SkeletonIDs = firstID;
while true
    nextID = -1;
    currentPoint = skeletonPoints(currentID).P;
    tempNeighborSkeletonIDs = skeletonPoints(currentID).NeighborSkeletonIDs;
    for i = 1:1:length(tempNeighborSkeletonIDs)
        tempNeighborSkeletonPoint = skeletonPoints(tempNeighborSkeletonIDs(i)).P;
        distance = sqrt(sum((currentPoint - tempNeighborSkeletonPoint).^2));
        if skeletonPoints(tempNeighborSkeletonIDs(i)).IsIgnore || (distance > parameters.distanceTracingThreshold)
            continue;
        end
        if distance <= parameters.nearThreshold
            skeletonPoints(tempNeighborSkeletonIDs(i)).P = [88888.8,88888.8,88888.8];
            skeletonPoints(tempNeighborSkeletonIDs(i)).IsFix = false;
            skeletonPoints(tempNeighborSkeletonIDs(i)).IsIgnore = true;
            skeletonPoints(tempNeighborSkeletonIDs(i)).IsBranch = false;
            skeletonPoints(tempNeighborSkeletonIDs(i)).NeighborSkeletonIDs = [];
            skeletonPoints(tempNeighborSkeletonIDs(i)).NeighborSkeletonPoints = [];
            skeletonPoints(tempNeighborSkeletonIDs(i)).NeighborTLSIDs = [];
            skeletonPoints(tempNeighborSkeletonIDs(i)).NeighborTLSPoints = [];
            skeletonPoints(tempNeighborSkeletonIDs(i)).NeighborTLSDensities = [];
            skeletonPoints(tempNeighborSkeletonIDs(i)).Sigma = 0;
            skeletonPoints(tempNeighborSkeletonIDs(i)).AverageTerm = [0,0,0];
            skeletonPoints(tempNeighborSkeletonIDs(i)).AverageWeightSum = 0;
            skeletonPoints(tempNeighborSkeletonIDs(i)).RepulsionTerm = [0,0,0];
            skeletonPoints(tempNeighborSkeletonIDs(i)).RepulsionWeightSum = 0;
            continue;
        end
        projDistance = (tempNeighborSkeletonPoint - currentPoint)*searchDirection';
        if projDistance >= 0
            nextID = tempNeighborSkeletonIDs(i);
            break;
        end
    end
    if nextID == -1
        break;
    end
    nextPoint = skeletonPoints(nextID).P;
    newDirection = (nextPoint - currentPoint)/sqrt(sum((nextPoint - currentPoint).^2));
    searchAngle = acos(newDirection*searchDirection')*180/pi;
    if (searchAngle > parameters.angleTracingThreshold) || ~skeletonPoints(nextID).IsFix || skeletonPoints(nextID).IsBranch
        break;
    else
        newBranch.SkeletonIDs = [newBranch.SkeletonIDs,nextID];
        searchDirection = newDirection;
        currentID = nextID;
    end
end

end