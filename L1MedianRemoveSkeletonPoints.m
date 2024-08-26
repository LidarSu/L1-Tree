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
% L1MedianRemoveSkeletonPoints.m    The function to ensure the spacing between skeleton points is not too close
% 
% Version 1.0
% Latest update     28 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
% 
% INPUTS:
% skeletonPoints    The struct that records the attributes of skeleton points
% parameters    The updated parameter set
% globalBar    The global waiting bar
% iterationNum    The number of iterations
% maxIterationNum    The maximum number of iterations
% 
% OUTPUTS:
% skeletonPoints    The updated struct that records the attributes of skeleton points
% ------------------------------------------------------------------------------

function skeletonPoints = L1MedianRemoveSkeletonPoints(skeletonPoints,parameters,globalBar,iterationNum,maxIterationNum)

waitbar(iterationNum/maxIterationNum,globalBar,'Remove too close skeleton points');

for i = 1:1:size(skeletonPoints,2)
    tempPoint = skeletonPoints(i).P;
    tempNeighborSkeletonIDs = skeletonPoints(i).NeighborSkeletonIDs;
    if isempty(tempNeighborSkeletonIDs)
        continue;
    end
    for j = 1:1:length(tempNeighborSkeletonIDs)
        if skeletonPoints(tempNeighborSkeletonIDs(j)).IsFix || skeletonPoints(tempNeighborSkeletonIDs(j)).IsIgnore
            continue;
        end
        tempNeighborSkeletonPoint = skeletonPoints(tempNeighborSkeletonIDs(j)).P;
        distance = sqrt(sum((tempNeighborSkeletonPoint - tempPoint).^2));
        if distance < parameters.nearThreshold
            skeletonPoints(tempNeighborSkeletonIDs(j)).P = [88888.8,88888.8,88888.8];
            skeletonPoints(tempNeighborSkeletonIDs(j)).IsFix = false;
            skeletonPoints(tempNeighborSkeletonIDs(j)).IsIgnore = true;
            skeletonPoints(tempNeighborSkeletonIDs(j)).IsBranch = false;
            skeletonPoints(tempNeighborSkeletonIDs(j)).NeighborSkeletonIDs = [];
            skeletonPoints(tempNeighborSkeletonIDs(j)).NeighborSkeletonPoints = [];
            skeletonPoints(tempNeighborSkeletonIDs(j)).NeighborTLSIDs = [];
            skeletonPoints(tempNeighborSkeletonIDs(j)).NeighborTLSPoints = [];
            skeletonPoints(tempNeighborSkeletonIDs(j)).NeighborTLSDensities = [];
            skeletonPoints(tempNeighborSkeletonIDs(j)).Sigma = 0;
            skeletonPoints(tempNeighborSkeletonIDs(j)).AverageTerm = [0,0,0];
            skeletonPoints(tempNeighborSkeletonIDs(j)).AverageWeightSum = 0;
            skeletonPoints(tempNeighborSkeletonIDs(j)).RepulsionTerm = [0,0,0];
            skeletonPoints(tempNeighborSkeletonIDs(j)).RepulsionWeightSum = 0;
        end
    end
end

end