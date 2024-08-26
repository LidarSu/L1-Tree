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
% L1TreeUpdate.m    The function for generating the struct of the 3D tree model
%
% Version 1.0
% Latest update     03 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% treeModel    The input 3D tree model
% skeletonPoints    Input skeleton points
% rootID    The ID of the root point
%
% OUTPUTS:
% treeModel    The updated 3D tree model
% skeletonPoints    Updated skeleton points
% ------------------------------------------------------------------------------

function [treeModel,skeletonPoints] = L1TreeUpdate(treeModel,skeletonPoints,rootID)

if ~isstruct(treeModel)
    treeGraph = graph(treeModel(:,1),treeModel(:,2),treeModel(:,3));
    [~,pred] = minspantree(treeGraph,'Method','dense','Type','tree','Root',rootID);
    rootID = find(pred == 0);
    nextID = find(pred == rootID);
    treeModel = [];
    currentCenCode = 'T';
    treeModel = L1TreeGetCenCode(treeModel,pred,currentCenCode,rootID,nextID);
end

%% Reorder skeleton points
occurrenceOrder = [];
for i = 1:1:size(treeModel,2)
    for j = 1:1:length(treeModel(i).SkeletonIDs)
        if ~ismember(treeModel(i).SkeletonIDs(j),occurrenceOrder)
            occurrenceOrder = [occurrenceOrder,treeModel(i).SkeletonIDs(j)];
        end
    end
end
for i = 1:1:size(treeModel,2)
    newSkeletonIDs = [];
    for j = 1:1:length(treeModel(i).SkeletonIDs)
        selectID = find(occurrenceOrder == treeModel(i).SkeletonIDs(j));
        newSkeletonIDs = [newSkeletonIDs,selectID];
    end
    treeModel(i).SkeletonIDs = newSkeletonIDs;
end
skeletonPoints = skeletonPoints(occurrenceOrder,:);

%% Output
updatedTreeModel = [];
for i = 1:1:size(treeModel,2)
    tempBranch.SkeletonIDs = treeModel(i).SkeletonIDs;
    tempBranch.Curve = skeletonPoints(treeModel(i).SkeletonIDs,:);
    tempBranch.CenCode = treeModel(i).CenCode;
    tempBranch.CenOrder = treeModel(i).CenOrder;
    tempBranch.CenRadius = 0;
    tempBranch.CenLength = 0;
    tempBranch.IsCenCredible = true;
    tempBranch.IsTerminal = true;
    for j = 1:1:size(treeModel,2)
        if contains(treeModel(j).CenCode,treeModel(i).CenCode) && (treeModel(j).CenOrder == treeModel(i).CenOrder+1)
            tempBranch.IsTerminal = false;
            break;
        end
    end
    tempBranch.GaCode = [];
    tempBranch.GaOrder = 0;
    tempBranch.GaRadius = 0;
    tempBranch.GaLength = 0;
    tempBranch.IsGaCredible = true;
    updatedTreeModel = [updatedTreeModel,tempBranch];
end
treeModel = updatedTreeModel;

end