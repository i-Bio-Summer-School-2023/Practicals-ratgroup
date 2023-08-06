function cv = crossvalPartition(n, kfold)
%Creates a partition of 1:n indices into k-fold.


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