function sparsityIndex = getSparsity(t, o)
% getSparsity Compute the sparsity index of a tuning curve.
%
% sparsityIndex = getSparsity(t, o) computes the selectivity of a tuning curve t
% given an occupancy o. The sparsity index is defined as 1 minus the squared sum
% of (t.*o)/sum(o) divided by the squared sum of (t.^2.*o)/sum(o), as proposed by
% Jung, Wiener, and McNaughton (JNS, 1994). A higher value indicates greater sparsity
% and selectivity in the tuning curve.
%
% INPUTS:
% - t: Tuning curve representing firing rates.
% - o: Occupancy map corresponding to the tuning curve.
%
% OUTPUT:
% - sparsityIndex: Sparsity index of the tuning curve.
%
% USAGE:
% sparsityIndex = getSparsity(t, o);
%
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%

%vectorizing to make sure t and o have the same shape
t = t(:);
o = o(:);

%Excluding values where either t or o are missing
valididx = ~isnan(t) & ~isnan(o);
t = t(valididx);
o = o(valididx);

%Calculating the sparsity index
sparsityIndex = 1 - sum(t(:) .* o(:) / sum(o(:))).^2 / sum(t(:).^2 .* o(:) / sum(o(:)));

end