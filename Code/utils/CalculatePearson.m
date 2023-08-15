function corr = CalculatePearson(t1,t2,d)
% CalculatePearson - Compute Pearson's correlation between Nd arrays t1 and t2.
%
%   corr = CalculatePearson(t1, t2, d) computes the Pearson's correlation
%   coefficient between the Nd arrays t1 and t2, integrating over dimensions
%   specified in d. If d is not provided, the function operates along the first
%   non-singleton dimension.
%
% Input:
%   t1, t2 - Input Nd arrays for correlation calculation.
%   d - Dimensions over which to compute the correlation (optional).
%
% Output:
%   corr - Pearson's correlation coefficient.
%
% Usage:
%   corr = CalculatePearson(tuning1, tuning2, 2);
%
% Note:
%   The arrays t1 and t2 should have the same dimensions.
%
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"

%%
sz = size(t1);
if nargin < 3
    d = find(sz > 1, 1, 'first');
end

if ~isequal(sz, size(t2))
    error('arrays should have the same dimension')
end

t1 = t1 - mean(t1,d,'omitnan');
t2 = t2 - mean(t2,d,'omitnan');
corr = sum(t1.*t2, d, 'omitnan')./sqrt(sum(t1.^2, d, 'omitnan').*sum(t2.^2, d, 'omitnan'));

end