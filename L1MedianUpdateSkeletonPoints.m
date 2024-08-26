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
% L1MedianUpdateSkeletonPoints.m    The function for calculating the average term of each skeleton point
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
% localBar    The local waiting bar
% 
% OUTPUTS:
% skeletonPoints    The updated struct that records the attributes of skeleton points
% averageMovingDistance    The average moving distance of skeleton points in each iteration
% ------------------------------------------------------------------------------

function [skeletonPoints,averageMovingDistance] = L1MedianUpdateSkeletonPoints(skeletonPoints,parameters,localBar)

waitbar(4/5,localBar,'Update skeleton points');

minSigma = min([skeletonPoints(:).Sigma]);
maxSigma = max([skeletonPoints(:).Sigma]);
rangeSigma = max(abs(maxSigma - minSigma),0.01);
minMu = parameters.minMu;
maxMu = parameters.maxMu;
rangeMu = abs(maxMu - minMu);

movingDistance = 0;
movingPointsNum = 0;
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsFix || skeletonPoints(i).IsIgnore
        continue;
    end
    tempPoint = skeletonPoints(i).P;
    tempSigma = skeletonPoints(i).Sigma;
    tempAverageWeightSum = skeletonPoints(i).AverageWeightSum;
    tempRepulsionWeightSum = skeletonPoints(i).RepulsionWeightSum;
    tempMu = (rangeMu/rangeSigma)*(tempSigma - minSigma) + minMu;
    if (tempAverageWeightSum > 0) && (tempRepulsionWeightSum > 0) && (tempMu >= 0)
        skeletonPoints(i).P = skeletonPoints(i).AverageTerm/tempAverageWeightSum + tempMu*skeletonPoints(i).RepulsionTerm/tempRepulsionWeightSum;
        movingDistance = movingDistance + sqrt(sum((skeletonPoints(i).P - tempPoint).^2));
        movingPointsNum = movingPointsNum + 1;
    end
end
averageMovingDistance = movingDistance/movingPointsNum;

end