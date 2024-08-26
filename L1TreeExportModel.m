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
% L1TreeExportModel.m    The function for export the 3D tree model
%
% Version 1.0
% Latest update     05 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% fileName     Filename for point cloud data without extensions
% treeModel    The 3D tree model
%
% OUTPUTS:
% outputFileName    The output filename
% ------------------------------------------------------------------------------

function outputFileName = L1TreeExportModel(fileName,treeModel)

maxCenOrder = max([treeModel(:).CenOrder]);
for i = 1:1:maxCenOrder
    for j = 1:1:size(treeModel,2)
        if (treeModel(j).CenOrder ~= i) || ~isnan(treeModel(j).CenRadius)
            continue;
        end
        for k = 1:1:size(treeModel,2)
            if contains(treeModel(j).CenCode,treeModel(k).CenCode) && (treeModel(k).CenOrder == treeModel(j).CenOrder-1)
                parentBranchID = k;
            end
        end
        treeModel(j).CenOrder = treeModel(parentBranchID).CenRadius*0.75;
    end
end

outputMatrix = [];
for i = 1:1:size(treeModel,2)
    tempBranch = treeModel(i);
    tempOutputMatrix = [ones(size(tempBranch.Curve,1),1)*i,tempBranch.Curve,ones(size(tempBranch.Curve,1),1)*tempBranch.CenOrder,...
                        ones(size(tempBranch.Curve,1),1)*tempBranch.CenRadius];
    outputMatrix = [outputMatrix;tempOutputMatrix];
end
outputFileName = ['./Results/',fileName,'_TreeModelDraw.txt'];
writematrix(outputMatrix,outputFileName,'Delimiter',' ');

end