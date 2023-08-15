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
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%
% Calculate Residual Sum of Squares (RSS)
RSS = sum((y - ypred).^2, 'omitnan');

% Calculate Mean of the original data (y)
m = mean(y, 'omitnan');

% Calculate Total Sum of Squares (TSS)
TSS = sum((y - m).^2, 'omitnan');

% Calculate Explained Variance (EV)
EV = 1 - (RSS/TSS);
end