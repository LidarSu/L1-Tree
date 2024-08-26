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
% L1MedianCalculateSigma.m    The function for calculating the directional degree of each skeleton point
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
% ------------------------------------------------------------------------------

function skeletonPoints = L1MedianCalculateSigma(skeletonPoints,parameters,localBar)

waitbar(1/5,localBar,'Calculate directional degree (sigma)');

for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsFix || skeletonPoints(i).IsIgnore
        continue;
    end
    tempPoint = skeletonPoints(i).P;
    tempNeighborSkeletonIDs = skeletonPoints(i).NeighborSkeletonIDs;
    if length(tempNeighborSkeletonIDs) <= 3
        skeletonPoints(i).Sigma = 0.5;
    else
        tempNeighborSkeletonPoints = skeletonPoints(i).NeighborSkeletonPoints;
        distances = sqrt(sum((tempNeighborSkeletonPoints - tempPoint).^2,2));
        weights = exp(distances.^2*(-1/(parameters.searchRange*2/2)^2));
        C = ((tempPoint - tempNeighborSkeletonPoints).*weights)'*(tempPoint - tempNeighborSkeletonPoints);
        [~,eigValues] = eig(C);
        eigValues = diag(eigValues);
        skeletonPoints(i).Sigma = max(eigValues)/sum(eigValues);
    end
end

end