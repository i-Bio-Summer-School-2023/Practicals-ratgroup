function [XMaxDec, YMaxDec, XMeanDec, YMeanDec] = ComputeBayesMAP(map, spk, dectimewin, mapprior)
% ComputeBayesMAP - Computes the decoded position using Bayesian Maximum
%                  A Posteriori (MAP) or Mean estimates.
%
%   [XMaxDec, YMaxDec, XMeanDec, YMeanDec] = ComputeBayesMAP(map, spk, dectimewin, mapprior)
%   computes the decoded position as the maximum or the mean of the
%   posterior probability distribution for 2D models. The posterior
%   probabilities are computed using a Bayesian approach, assuming
%   independence between spike trains across cells. The function requires
%   the tuning curve maps, spike counts, decoding time window, and optional
%   prior distribution.
%
% INPUTS:
% - map: 2D tuning curves (ncells x nYbins x nXbins in spikes/second).
% - spk: Spike counts across time for all cells (ntimes x ncells).
% - dectimewin: Decoding window in seconds, corresponding to the time
%   window used for spike counts.
% - mapprior (optional): Prior distribution of decoded variables. Default
%                        is a flat prior.
%
% OUTPUTS:
% - XMaxDec: Decoded X position using the maximum of the posterior.
% - YMaxDec: Decoded Y position using the maximum of the posterior.
% - XMeanDec: Decoded X position using the mean of the posterior.
% - YMeanDec: Decoded Y position using the mean of the posterior.
%
% USAGE:
% [XMaxDec, YMaxDec, XMeanDec, YMeanDec] = ComputeBayesMAP(map, spk, dectimewin, mapprior);
%
% Written by J Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%

if nargin < 4
    %if mapprior is not provided, we'll use a flat prior. Units don't
    %matter as this will be normalized later
    mapprior = ones([1 size(map, [2 3])]);
end

%Permuting dimensions to get first dimension as a singleton and do the 
%element-wise multiplication with spk more conveniently
map = permute(map, [4 1 2 3]);
mapprior = permute(mapprior, [4 1 2 3]);

%Computing the posterior probability P(X,Y | spk). We'll do it as a for
%loop over the number of cells to avoid memory issues.
ncells = size(spk, 2);
map = map + eps;%to avoid reaching precision limit
Posterior = mapprior .* exp(-dectimewin .* sum(map, 2));
for icell = 1:ncells
    Posterior = Posterior .* (map(1,icell,:,:) .^ spk(:,icell));
end

%We should end up with a matrix of probability of size Ntimes x nXbins x
%nYbins
Posterior = permute(Posterior, [1 3 4 2]);

%Normalizing so to that the sum of probabilities over positions equals 1.
Posterior = Posterior ./ sum(Posterior, [2 3], 'omitnan');

%Taking the decoded position as the maximum of the posterior probability
%distribution (M.A.P. estimate)
[~, DecMax] = max(Posterior, [], [2 3], 'linear');

%Converting DecMax from linear indices to subscripts
[~,YMaxDec, XMaxDec] = ind2sub(size(Posterior), DecMax);

%Taking the decoded position as the expected value of the position given 
%its posterior probability distribution. We'll compute the mean on the
%marginal of the posterior distribution along the considered variable.
%First for the Y variable.
YPosterior = squeeze(sum(Posterior, 3, 'omitnan')) / 2;
Ybins = (1:size(Posterior, 2))';
YPosterior(isnan(YPosterior)) = 0;
YMeanDec = (YPosterior * Ybins) ./ sum(YPosterior, 2);

%Then for the X variable.
XPosterior = squeeze(sum(Posterior, 2, 'omitnan')) / 2;
Xbins = (1:size(Posterior, 3))';
XPosterior(isnan(XPosterior)) = 0;
XMeanDec = (XPosterior * Xbins) ./ sum(XPosterior, 2);

%Ignoring decoded positions if no cell fired (optional)
mua = sum(spk, 2);
XMaxDec(mua == 0) = NaN;
YMaxDec(mua == 0) = NaN;
XMeanDec(mua == 0) = NaN;
YMeanDec(mua == 0) = NaN;
end