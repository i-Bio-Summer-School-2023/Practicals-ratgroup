function [SInfoperspike, SInfo] = getSpatialinfo(t, o)
% getSpatialinfo computes the spatial information in a tuning curve.
%
% SInfoperspike = getSpatialinfo(t, o) computes the spatial information in a tuning curve t,
% considering an occupancy o. The result is returned in bits per spike.
%
% INPUTS:
% - t: Tuning curve (often representing firing rates).
% - o: Occupancy map corresponding to the tuning curve.
%
% OUTPUT:
% - SInfoperspike: Spatial information in bits per spike.
% - SInfo: Spatial information in bits per second.
%
% USAGE:
% SInfoperspike = getSpatialinfo(t, o);
%
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%
%%
%vectorizing to make sure t and o have the same shape
t = t(:);
o = o(:);

%Excluding values where either t or o are missing
valididx = ~isnan(t) & ~isnan(o);
t = t(valididx);
o = o(valididx);

%Mean rate of the cell
meanRate = sum(t(:) .* o(:) / sum(o(:)));

%Spatial information in bits per seconds.
SInfo = sum(o(:) / sum(o(:)) .* t(:) .* log2(t(:) / meanRate));

%Converting in bits per spike.
SInfoperspike = SInfo/meanRate;
end