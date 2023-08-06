function [LLH, BIC, AIC] = computeLLH_poisson(y, ypred, k)
%Compute the log likelihood for a Poisson model. y is the original signal;
%ypred, the model prediction and k is the total number of model parameters.
%k is only necessary if the Bayesian Information Criterion and Akaike
%Information Criterion are required.

%probability density function for a normal distribution with std = s.
pd = poisspdf(y, ypred);
pd(pd == 0) = eps;

%Log likelihood is the sum of the log of the probabilities
LLH = nansum(log(pd));

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