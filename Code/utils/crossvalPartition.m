function cv = crossvalPartition(n, kfold)
% cv = crossvalPartition(n, kfold) returns a structure containing k-fold cross-validation
% partition sets where training and test sets are contiguous within each partition.
%
% INPUTS:
% - n: Total number of data points.
% - kfold: Number of folds for cross-validation.
%
% OUTPUT:
% - cv: Struct with fields trainsets and testsets, each containing kfold cell arrays
%   defining the training and testing indices for each fold.
%
% USAGE:
% cv = crossvalPartition(n, kfold);
%
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%
cv.trainsets = cell(1,kfold);
cv.testsets = cell(1,kfold);

%Edges of the subparts, taken as contiguous.
kidx = floor(n./kfold).*(0:kfold);
kidx(end) = n;

%Defining trainsets and tests as boolean vectors within the range of the
%partition edges.
for k = 1:kfold
    cv.trainsets{k} = false(n,1);
    cv.testsets{k} = false(n,1);
    
    cv.testsets{k}((kidx(k)+1):kidx(k+1)) = true;
    cv.trainsets{k}(~cv.testsets{k}) = true;
end
end