function [LLH, BIC, AIC] = computeLLH_normal(y, ypred, k)
% computeLLH_normal Compute log likelihood, Bayesian Information Criterion (BIC),
% and Akaike Information Criterion (AIC) for a Gaussian model.
%
% [LLH, BIC, AIC] = computeLLH_normal(y, ypred, k) computes the log likelihood (LLH)
% for a Gaussian model given the original signal (y) and its model prediction (ypred).
% The total number of model parameters (k) is optionally provided for calculating BIC and AIC.
%
% INPUTS:
% - y: Original signal.
% - ypred: Model prediction.
% - k: Total number of model parameters (optional for BIC and AIC).
%
% OUTPUTS:
% - LLH: Log likelihood of the Gaussian model.
% - BIC: Bayesian Information Criterion.
% - AIC: Akaike Information Criterion.
%
% USAGE:
% [LLH, BIC, AIC] = computeLLH_normal(y, ypred, k);
%
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%
%std of the model
s = nanstd((y - ypred), 1);

%probability density function for a normal distribution with std = s.
pdfun = normpdf(y, ypred, s);
pdfun(pdfun == 0) = eps;

%Log likelihood is the sum of the log of the probabilities
LLH = nansum(log(pdfun));

%Alternatively, one could compute the log likelihood by hand as follows
% n = numel(y);
% LLH = -n/2 * log(2*pi) -n/2 * log(s^2) -1/(2*s^2)*sum((y - ypred).^2);

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