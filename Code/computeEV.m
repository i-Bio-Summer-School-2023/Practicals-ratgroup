function EV = computeEV(y, ypred)
% computeEV computes the Explained Variance (EV) of a model's prediction.
%
% EV = computeEV(y, ypred) computes the Explained Variance (EV) of a model's prediction
% given the original data (y) and its prediction (ypred).
%
% INPUTS:
% - y: Original data.
% - ypred: Model prediction.
%
% OUTPUT:
% - EV: Explained Variance.
%
% USAGE:
% EV = computeEV(y, ypred);
%
%
% Written by J.Fournier 08/2023 for the iBio Summer school


RSS = sum((y - ypred).^2, 'omitnan');
m = mean(y, 'omitnan');
TSS = sum((y - m).^2, 'omitnan');

EV = 1 - (RSS/TSS);
end