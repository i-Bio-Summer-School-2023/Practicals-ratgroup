function pred = getmapPred(map, X, scaling)
% return the output of a nD model (with n < or  = to 3 for now...)

lintuning = map(:);
Xinput = [];
for k = 1:size(X,2)
    Xinput = cat(2, Xinput, {X(:,k)});
end
linidx = sub2ind(size(map), Xinput{:});
scaling(scaling == 0) = NaN;
pred = NaN(size(linidx));
pred(~isnan(linidx))  = lintuning(linidx(~isnan(linidx))) .* scaling(~isnan(linidx));
end