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
% Written by J.Fournier 08/2023 for the iBio Summer school


sparsityIndex = 1 - nansum(t(:).*o(:)/nansum(o(:))).^2/nansum(t(:).^2.*o(:)/nansum(o(:)));

end