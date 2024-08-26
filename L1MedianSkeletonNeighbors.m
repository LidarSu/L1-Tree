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
% L1MedianSkeletonNeighbors.m    The function for identifying neighbors of each skeleton point
% 
% Version 1.0
% Latest update     28 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
% 
% INPUTS:
% skeletonPoints    The struct that records the attributes of skeleton points
% tlsPoints    The struct that records the attributes of TLS points
% parameters    The updated parameter set
% localBar    The local waiting bar
% 
% OUTPUTS:
% skeletonPoints    The updated struct that records the attributes of skeleton points
% ------------------------------------------------------------------------------

function skeletonPoints = L1MedianSkeletonNeighbors(skeletonPoints,tlsPoints,parameters,localBar)

waitbar(0,localBar,'Identify neighbors of skeleton points');

%% Load data
skeletonPointPs = [skeletonPoints(:).P];
skeletonPointPs = reshape(skeletonPointPs,3,size(skeletonPoints,2))';
tlsPointPs = [tlsPoints(:).P];
tlsPointPs = reshape(tlsPointPs,3,size(tlsPoints,2))';
densities = [tlsPoints(:).Density];
isIgnore = [skeletonPoints(:).IsIgnore];
matchTable = find(~isIgnore);

%% Identify skeleton neighbors of each skeleton point
[skeletonNeighborIDs_1,skeletonDistances_1] = rangesearch(skeletonPointPs(~isIgnore,:),skeletonPointPs(~isIgnore,:),parameters.baseRadius*2,'Distance','euclidean','NSMethod','kdtree');
skeletonNeighborIDs_2 = rangesearch(skeletonPointPs(~isIgnore,:),skeletonPointPs(~isIgnore,:),parameters.searchRange*2,'Distance','euclidean','NSMethod','kdtree');
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsIgnore
        skeletonPoints(i).NeighborSkeletonPoints = [];
        continue;
    end
    if skeletonPoints(i).Location/skeletonPoints(i).PathLength*parameters.baseRadius*2 <= parameters.searchRange*2
        selectID = find(matchTable == i);
        skeletonPoints(i).NeighborSkeletonIDs = setdiff(skeletonNeighborIDs_2{selectID},selectID);
        skeletonPoints(i).NeighborSkeletonIDs = matchTable(skeletonPoints(i).NeighborSkeletonIDs);
        skeletonPoints(i).NeighborSkeletonPoints = skeletonPointPs(skeletonPoints(i).NeighborSkeletonIDs,:);
    else
        selectID = find(matchTable == i);
        neighborIDs = skeletonNeighborIDs_1{selectID};
        neighborDistances = skeletonDistances_1{selectID};
        selectIDs = find(neighborDistances <= skeletonPoints(i).Location/skeletonPoints(i).PathLength*parameters.baseRadius*2);
        skeletonPoints(i).NeighborSkeletonIDs = setdiff(neighborIDs(selectIDs),selectID);
        skeletonPoints(i).NeighborSkeletonIDs = matchTable(skeletonPoints(i).NeighborSkeletonIDs);
        skeletonPoints(i).NeighborSkeletonPoints = skeletonPointPs(skeletonPoints(i).NeighborSkeletonIDs,:);
    end
end

%% Identify TLS neighbors if each skeleton point
[tlsNeighborsIDs_1,tlsDistances_1] = rangesearch(tlsPointPs,skeletonPointPs(~isIgnore,:),parameters.baseRadius*2,'Distance','euclidean','NSMethod','kdtree');
tlsNeighborsIDs_2 = rangesearch(tlsPointPs,skeletonPointPs(~isIgnore,:),parameters.searchRange*2,'Distance','euclidean','NSMethod','kdtree');
for i = 1:1:size(skeletonPoints,2)
    if skeletonPoints(i).IsIgnore
        skeletonPoints(i).NeighborTLSPoints = [];
        skeletonPoints(i).NeighborTLSDensities = [];
        continue;
    end
    if skeletonPoints(i).Location/skeletonPoints(i).PathLength*parameters.baseRadius*2 <= parameters.searchRange*2
        selectID = find(matchTable == i);
        skeletonPoints(i).NeighborTLSIDs = tlsNeighborsIDs_2{selectID};
        skeletonPoints(i).NeighborTLSPoints = tlsPointPs(skeletonPoints(i).NeighborTLSIDs,:);
        skeletonPoints(i).NeighborTLSDensities = densities(skeletonPoints(i).NeighborTLSIDs);
    else
        selectID = find(matchTable == i);
        neighborIDs = tlsNeighborsIDs_1{selectID};
        neighborDistances = tlsDistances_1{selectID};
        selectIDs = find(neighborDistances <= skeletonPoints(i).Location/skeletonPoints(i).PathLength*parameters.baseRadius*2);
        skeletonPoints(i).NeighborTLSIDs = neighborIDs(selectIDs);
        skeletonPoints(i).NeighborTLSPoints = tlsPointPs(skeletonPoints(i).NeighborTLSIDs,:);
        skeletonPoints(i).NeighborTLSDensities = densities(skeletonPoints(i).NeighborTLSIDs);
    end
end

end