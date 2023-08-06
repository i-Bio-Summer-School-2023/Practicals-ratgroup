function [LLH, BIC, AIC] = computeLLH_normal(y, ypred, k)
%Compute the log likelihood for a Gaussian model. y is the original signal;
%ypred, the model prediction and k is the total number of model parameters.
%k is only necessary if the Bayesian Information Criterion and Akaike
%Information Criterion are required.

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