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
% L1MedianExtendBranchFromTail.m    Extend the branch from the end
% 
% Version 1.0
% Latest update     28 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
% 
% INPUTS:
% newBranch    The input branch which includes the IDs of its constituent skeleton points 
% skeletonPoints    The struct that records the attributes of skeleton points
% parameters    The updated parameter set
% 
% OUTPUTS:
% newBranch    The updated branch
% skeletonPoints    The updated struct that records the attributes of skeleton points
% ------------------------------------------------------------------------------

function [newBranch,skeletonPoints] = L1MedianExtendBranchFromTail(newBranch,skeletonPoints,parameters)

while true
    % Identify the end direction
    endPoint = skeletonPoints(newBranch.SkeletonIDs(end)).P;
    predPoint = skeletonPoints(newBranch.SkeletonIDs(end-1)).P;
    endDirection = (endPoint - predPoint)/sqrt(sum((endPoint - predPoint).^2));
    tempNeighborSkeletonIDs = skeletonPoints(newBranch.SkeletonIDs(end)).NeighborSkeletonIDs;
    % Identify the next point
    nextID = -1;
    for i = 1:1:length(tempNeighborSkeletonIDs)
        tempNeighborSkeletonPoint = skeletonPoints(tempNeighborSkeletonIDs(i)).P;
        distance = sqrt(sum((tempNeighborSkeletonPoint - endPoint).^2));
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
        projDistance = (tempNeighborSkeletonPoint - endPoint)*endDirection';
        if projDistance >= 0
            nextID = tempNeighborSkeletonIDs(i);
            nextSigma = skeletonPoints(tempNeighborSkeletonIDs(i)).Sigma;
            break;
        end
    end
    if (nextID == -1) || (nextSigma < parameters.sigmaTracingThreshold)
        break;
    end
    nextPoint = skeletonPoints(nextID).P;
    predPredPoint = skeletonPoints(newBranch.SkeletonIDs(end-2)).P;
    vector_1 = (predPoint - predPredPoint)/sqrt(sum((predPoint - predPredPoint).^2));
    vector_2 = (endPoint - predPoint)/sqrt(sum((endPoint - predPoint).^2));
    vector_3 = (nextPoint - endPoint)/sqrt(sum((nextPoint - endPoint).^2));
    angle_1 = acos(vector_1*vector_2')*180/pi;
    angle_2 = acos(vector_2*vector_3')*180/pi;
    if (angle_1 > parameters.angleTracingThreshold) || (angle_2 > parameters.angleTracingThreshold)
        break;
    end
    if skeletonPoints(nextID).IsBranch
        break;
    end
    % Update skeletonPoints.IsFix, skeletonPoints.IsBranch, and newBranch
    skeletonPoints(newBranch.SkeletonIDs(end)).IsFix = true;
    skeletonPoints(newBranch.SkeletonIDs(end)).IsBranch = true;
    skeletonPoints(nextID).IsFix = false;
    skeletonPoints(nextID).IsBranch = false;
    newBranch.SkeletonIDs = [newBranch.SkeletonIDs,nextID];
end

end