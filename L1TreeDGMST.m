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
% L1TreeDGMST.m    The function for extracting the minimum spanning tree from a directed graph
%                  Please see "https://wendy-xiao.github.io/posts/2020-07-10-chuliuemdond_algorithm/" for the principle of this algorithm
%
% Version 1.0
% Latest update     18 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% E    The list of edges
% rootID    The ID of the root point
% 
% OUTPUTS:
% E    The constituent edge of the minimum spanning tree
% ------------------------------------------------------------------------------

function E = L1TreeDGMST(E,rootID)

% Step 1: Reverse the graph and remove all edges starting from the root
rE = E;
rE(:,1) = E(:,2);
rE(:,2) = E(:,1);
delIDs = find(rE(:,1) == rootID);
rE(delIDs,:) = [];

% Step 2: Generate the subgraph with maximum outgoing edges
mE = [];
endPointIDs = unique(rE(:,1));
for i = 1:1:length(endPointIDs)
    endPointID = endPointIDs(i);
    selectIDs = find(rE(:,1) == endPointID);
    tempE = rE(selectIDs,:);
    [maxValue,maxID] = max(tempE(:,3));
    startPointID = tempE(maxID,2);
    mE = [mE;[endPointID,startPointID,maxValue]];
end

% Step 3: Determine whether there is a loop in the subgraph
startPointIDs = unique(mE(:,1));
for i = 1:1:length(startPointIDs)
    startPointID = startPointIDs(i);
    isVisited = [];
    stack = startPointID;
    isGetCycle = false;
    cyclePointIDs = [];
    while ~isempty(stack)
        currentPointID = stack(1);
        stack(1) = [];
        if ismember(currentPointID,isVisited)
            cyclePointIDs = [];
            while ~ismember(currentPointID,cyclePointIDs)
                cyclePointIDs = [cyclePointIDs,currentPointID];
                selectIDs = find(mE(:,1) == currentPointID);
                currentPointID = mE(selectIDs(1),2);
            end
            isGetCycle = true;
        end
        isVisited = [isVisited,currentPointID];
        if ismember(currentPointID,startPointIDs)
            selectIDs = find(mE(:,1) == currentPointID);
            stack = [stack,mE(selectIDs,2)'];
        end
        if isGetCycle
            break;
        end
    end
    if isGetCycle
        break;
    end
end

% Step 4: If there is no loop, invert the graph and return
if isempty(cyclePointIDs)
    E = mE;
    E(:,1) = mE(:,2);
    E(:,2) = mE(:,1);
    return
end

% Step 5: If there is a loop, break it
shrinkPointID = max(unique([E(:,1)',E(:,2)'])) + 1;
nE = [];
for i = 1:1:size(E,1)
    startPointID = E(i,1);
    endPointID = E(i,2);
    if ~ismember(startPointID,cyclePointIDs) && ismember(endPointID,cyclePointIDs)
        selectID = find(mE(:,1) == endPointID);
        newWeight = E(i,3) - mE(selectID,3);
        nE = [nE;[startPointID,shrinkPointID,newWeight,endPointID,E(i,3)]];
    end
    if ismember(startPointID,cyclePointIDs) && ~ismember(endPointID,cyclePointIDs)
        nE = [nE;[shrinkPointID,endPointID,E(i,3),startPointID,E(i,3)]];
    end
    if ~ismember(startPointID,cyclePointIDs) && ~ismember(endPointID,cyclePointIDs)
        nE = [nE;[startPointID,endPointID,E(i,3),-1,E(i,3)]];
    end
end
uniqueNE = unique(nE(:,1:2),'rows');
newWeights = zeros(size(uniqueNE,1),1);
oriPointIDs = zeros(size(uniqueNE,1),1);
oldWeights = zeros(size(uniqueNE,1),1);
for i = 1:1:size(uniqueNE,1)
    selectIDs = find((nE(:,1) == uniqueNE(i,1)) & (nE(:,2) == uniqueNE(i,2)));
    tempE = nE(selectIDs,:);
    [~,maxID] = max(tempE(:,3));
    newWeights(i) = tempE(maxID,3);
    oriPointIDs(i) = tempE(maxID,4);
    oldWeights(i) = tempE(maxID,5);
end
uniqueNE = [uniqueNE,newWeights];

% Step 6: Recursion
E = L1TreeDGMST(uniqueNE,rootID);
for i = 1:1:size(E,1)
    if E(i,1) == shrinkPointID
        startPointID = shrinkPointID;
        endPointID = E(i,2);
        selectID = find((uniqueNE(:,1) == startPointID) & (uniqueNE(:,2) == endPointID));
        E(i,1) = oriPointIDs(selectID);
        E(i,3) = oldWeights(selectID);
    end
    if E(i,2) == shrinkPointID
        startPointID = E(i,1);
        endPointID = shrinkPointID;
        selectID = find((uniqueNE(:,1) == startPointID) & (uniqueNE(:,2) == endPointID));
        E(i,2) = oriPointIDs(selectID);
        E(i,3) = oldWeights(selectID);
    end
end
addE = [];
for i = 1:1:size(mE,1)
    endPointID = mE(i,1);
    startPointID = mE(i,2);
    if ismember(startPointID,cyclePointIDs) && ismember(endPointID,cyclePointIDs)
        selectID = find(E(:,2) == endPointID);
        if isempty(selectID)
            addE = [addE;[startPointID,endPointID,mE(i,3)]];
        end
    end
end
E = [E;addE];