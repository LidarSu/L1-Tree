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
% L1TreeMovingTtest.m    The function for pruning terminal branches
%
% Version 1.0
% Latest update     04 Feb 2024
% 
% Copyright (C) 2024 Sulab, Institude of Botany, The Chinese Academy of Sciences
% If you have any questions about using the code, please contact Yuhao Feng (fengyuhao@pku.edu.cn)
%
% INPUTS:
% sequence    The input radius sequence
% 
% OUTPUT:
% breakID    The location of the break point
% ------------------------------------------------------------------------------

function breakID = L1TreeMovingTtest(sequence)

windowSize = max(floor(length(sequence)/10),5);

pValues = zeros(1,(length(sequence)-2*windowSize+1));
for i = windowSize:1:(length(sequence)-windowSize)
    sequence_1 = sequence((i-windowSize+1):i);
    sequence_1 = sequence_1(~isnan(sequence_1));
    sequence_2 = sequence((i+1):(i+windowSize));
    sequence_2 = sequence_2(~isnan(sequence_2));
    mean_1 = mean(sequence_1);
    mean_2 = mean(sequence_2);
    sd_1 = std(sequence_1);
    sd_2 = std(sequence_2);
    totalSD = sqrt((length(sequence_1)*sd_1^2 + length(sequence_2)*sd_2^2)/(length(sequence_1) + length(sequence_2)));
    tStatistics = (mean_2 - mean_1)/(totalSD*sqrt(1/length(sequence_1) + 1/length(sequence_2)));
    pValues(i-windowSize+1) = 1-tcdf(tStatistics,(length(sequence_1)+length(sequence_2)-2));
end

selectIDs = find(pValues < 0.01);
if isempty(selectIDs)
    breakID = [];
else
    selectIDs = selectIDs + windowSize - 1;
    breakID = selectIDs(1);
end

end