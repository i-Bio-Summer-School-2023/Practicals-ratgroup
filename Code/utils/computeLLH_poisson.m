function [LLH, BIC, AIC] = computeLLH_poisson(y, ypred, k)
% computeLLH_poisson Compute log likelihood, Bayesian Information Criterion (BIC),
% and Akaike Information Criterion (AIC) for a Poisson model.
%
% [LLH, BIC, AIC] = computeLLH_poisson(y, ypred, k) computes the log likelihood (LLH)
% for a Poisson model given the original signal (y) and its model prediction (ypred).
% The total number of model parameters (k) is optionally provided for calculating BIC and AIC.
%
% INPUTS:
% - y: Original signal.
% - ypred: Model prediction.
% - k: Total number of model parameters (optional for BIC and AIC).
%
% OUTPUTS:
% - LLH: Log likelihood of the Poisson model.
% - BIC: Bayesian Information Criterion.
% - AIC: Akaike Information Criterion.
%
% USAGE:
% [LLH, BIC, AIC] = computeLLH_poisson(y, ypred, k);
%
% See also: poisspdf
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%
%probability density function for a poisson distribution.
pd = poisspdf(y, ypred);
pd(pd == 0) = eps;

%Log likelihood is the sum of the log of the probabilities
LLH = sum(log(pd), 'omitnan');

%if k is provided, we also compute BIC and AIC of the model.
if nargin > 2
    N = numel(y);
    BIC = k * log(N) - 2 * LLH;
    AIC = 2 * k - 2 * LLH;
else
    BIC = NaN;
    AIC = NaN;
end
end