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
% L1MedianSmoothSigma.m    The function for smoothing the directional degree of each skeleton point
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

function skeletonPoints = L1MedianSmoothSigma(skeletonPoints,parameters,globalBar,iterationNum,maxIterationNum)

waitbar(iterationNum/maxIterationNum,globalBar,'Smooth directional degree (sigma)');

%% Update skeletonPoints.NeighborSkeletonIDs and skeletonPoints.NeighborSkeletonPoints
skeletonPointPs = [skeletonPoints(:).P];
skeletonPointPs = reshape(skeletonPointPs,3,size(skeletonPoints,2))';
[~,kIDs] = pdist2(skeletonPointPs,skeletonPointPs,'euc','Smallest',parameters.kSmoothSigma+1);
kIDs = kIDs(2:end,:);
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsFix || skeletonPoints(i).IsIgnore
        continue;
    end
    skeletonPoints(i).NeighborSkeletonIDs = [skeletonPoints(kIDs(:,i)).ID];
    tempNeighborSkeletonPoints = [skeletonPoints(kIDs(:,i)).P];
    skeletonPoints(i).NeighborSkeletonPoints = reshape(tempNeighborSkeletonPoints,3,length(skeletonPoints(i).NeighborSkeletonIDs))';
end

%% Recalculate and smooth the directional degree of each skeleton point
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsFix || skeletonPoints(i).IsIgnore
        continue;
    end
    tempNeighborSkeletonIDs = skeletonPoints(i).NeighborSkeletonIDs;
    if length(tempNeighborSkeletonIDs) <= 3
        skeletonPoints(i).Sigma = 0.5;
        continue;
    end
    selectIDs = [];
    for j = 1:1:length(tempNeighborSkeletonIDs)
        if ~skeletonPoints(tempNeighborSkeletonIDs(j)).IsFix && ~skeletonPoints(tempNeighborSkeletonIDs(j)).IsIgnore
            selectIDs = [selectIDs,j];
        end
    end
    skeletonPoints(i).NeighborSkeletonIDs = skeletonPoints(i).NeighborSkeletonIDs(selectIDs);
    skeletonPoints(i).NeighborSkeletonPoints = skeletonPoints(i).NeighborSkeletonPoints(selectIDs,:);
    tempPoint = skeletonPoints(i).P;
    tempNeighborSkeletonIDs = skeletonPoints(i).NeighborSkeletonIDs;
    tempNeighborSkeletonPoints = skeletonPoints(i).NeighborSkeletonPoints;
    if length(tempNeighborSkeletonIDs) <= 3
        skeletonPoints(i).Sigma = 0.95;
        continue;
    end
    C = (tempNeighborSkeletonPoints - tempPoint)'*(tempNeighborSkeletonPoints - tempPoint);
    [~,eigValues] = eig(C);
    eigValues = diag(eigValues);
    skeletonPoints(i).Sigma = max(eigValues)/sum(eigValues);
end

%% Smooth sigma
sigma = [skeletonPoints(:).Sigma];
sigma = 1 - sigma;
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsFix || skeletonPoints(i).IsIgnore
        continue;
    end
    tempSigma = sigma(i);
    tempNeighborSkeletonIDs = skeletonPoints(i).NeighborSkeletonIDs;
    tempSigmas = sigma(tempNeighborSkeletonIDs);
    skeletonPoints(i).Sigma = 1 - (tempSigma + sum(tempSigmas))/(1 + length(tempNeighborSkeletonIDs));
    if skeletonPoints(i).Sigma < 0
        skeletonPoints(i).Sigma = 0.5;
    end
end

%% Identify fixed skeleton points
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsBranch || skeletonPoints(i).IsIgnore
        continue;
    end
    if skeletonPoints(i).Sigma > parameters.sigmaThreshold
        skeletonPoints(i).IsFix = true;
    else
        skeletonPoints(i).IsFix = false;
    end
end

end