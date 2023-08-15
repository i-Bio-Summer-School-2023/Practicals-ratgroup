function [m, sd, r] = ComputeTriggeredAverage(R, S, idxwin, w)
% [m, sd, r] = ComputeTriggeredAverage(R, S, idxwin, w)
%
%Computes a triggered average of vector R based on timestamps in S over a 
%window of indices idxwin, weighted by w.
%
% INPUTS:
% - R: 1D array from which to compute the average.
% - S: list of indices from which to extract value of R.
% - idxwin: list of indices around values of S.
% - w (optional): list of the same size as S to weight R's snippets before
% averaging.
%
% OUTPUTS:
% - m: average of R triggered on indices in S on a window defined by
% idxwin.
% - sd: standard deviation of the average m.
% - r: snippets of R triggered on S. Each line correspond to one snippet.
%
% USAGE:
% [m, sd, r] = ComputeTriggeredAverage(R, S, idxwin, [w]);
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%

if nargin < 4
    w = ones(size(S));
end

w = w(:);

R = R(:);

%Padding R
idxmax = max(abs(idxwin));
R = [NaN(idxmax,1) ; R; NaN(idxmax,1)];

%Extracting snippets of R around timestamps in S
r = arrayfun(@(x) R(x + idxwin), S + idxmax, 'UniformOutput', false);
r = cat(2, r{:})';
r = r .* w;

%Average across snippets
m = mean(r, 1, 'omitnan');

%s.d. across snippets
sd = std(r, 0, 1, 'omitnan');

end