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
% L1MedianCalculateDensity.m    The function for calculating the point density of each TLS point
% 
% Version 1.0
% Latest update     28 Jan 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% tlsPoints    The struct that records the attributes of TLS points
% parameters    The updated parameter set
% globalBar    The global waiting bar
% iterationNum    The number of iterations
% maxIterationNum    The maximum number of iterations
% 
% OUTPUTS:
% tlsPoints    The updated struct that records the attributes of TLS points
% ------------------------------------------------------------------------------

function tlsPoints = L1MedianCalculateDensity(tlsPoints,parameters,globalBar,iterationNum,maxIterationNum)

waitbar(iterationNum/maxIterationNum,globalBar,'Identify neighbors of TLS points');

tlsPointPs = [tlsPoints(:).P];
tlsPointPs = reshape(tlsPointPs,3,size(tlsPoints,2))';
[tlsPointIDs,distances] = rangesearch(tlsPointPs,tlsPointPs,parameters.initialSearchRange*2,'Distance','euclidean','NSMethod','kdtree');
for i = 1:size(tlsPoints,2)
    tlsPoints(i).NeighborTLSIDs = setdiff(tlsPointIDs{i},i);
    tlsPoints(i).NeighborTLSPoints = tlsPointPs(tlsPoints(i).NeighborTLSIDs,:);
    if isempty(tlsPoints(i).NeighborTLSIDs)
        tlsPoints(i).Density = 1;
    else
        tempDistances = setdiff(distances{i},0);
        weights = exp(tempDistances.^2*(-1/(parameters.initialSearchRange*2/2)^2));
        tlsPoints(i).Density = 1/(1 + sum(weights));
    end
end

end