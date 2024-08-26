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
% L1TreeGetCenCode.m    The function for identifying the centrifugal code of each branch
%
% Version 1.0
% Latest update     02 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% treeModel    The tree model where the centrifugal order of each branch has been identified
% pred    The list of predecessor IDs for each skeleton point
% currentCenCode    The centrifugal code of the current branch
% rootID    The ID of the root point
% nextID    The ID of the point next to the root point
%
% OUTPUTS:
% treeModel    The updated tree model
% ------------------------------------------------------------------------------

function treeModel = L1TreeGetCenCode(treeModel,pred,currentCenCode,rootID,nextID)

% Recursive calculation
newBranch.SkeletonIDs = [rootID,nextID];
newBranch.CenCode = currentCenCode;
newBranch.CenOrder = length(strfind(currentCenCode,'_')) + 1;
while true
    currentID = newBranch.SkeletonIDs(end);
    nextIDs = find(pred == currentID);
    if isempty(nextIDs)
        treeModel = [treeModel,newBranch];
        break;
    elseif length(nextIDs) == 1
        nextID = nextIDs(1);
        newBranch.SkeletonIDs = [newBranch.SkeletonIDs,nextID];
    else
        treeModel = [treeModel,newBranch];
        for i = 1:1:length(nextIDs)
            rootID = currentID;
            nextID = nextIDs(i);
            nextCenCode = [currentCenCode,'_B',num2str(i)];
            treeModel = L1TreeGetCenCode(treeModel,pred,nextCenCode,rootID,nextID);
        end
        break;
    end
end

end