function SInfoperspike = getSpatialInfo(t, o)
% getSpatialInfo computes the spatial information in a tuning curve.
%
% SInfoperspike = getSpatialInfo(t, o) computes the spatial information in a tuning curve t,
% considering an occupancy o. The result is returned in bits per spike.
%
% INPUTS:
% - t: Tuning curve (often representing firing rates).
% - o: Occupancy map corresponding to the tuning curve.
%
% OUTPUT:
% - SInfoperspike: Spatial information in bits per spike.
%
% USAGE:
% SInfoperspike = getSpatialInfo(t, o);
%
%
% Written by J.Fournier 08/2023 for the iBio Summer school


%Mean rate of the cell
meanRate = nansum(t(:) .* o(:) / nansum(o(:)));

%Spatial information in bits per seconds.
SInfo = nansum(o(:) / nansum(o(:)) .* t(:) .* log2(t(:) / meanRate));

%Converting in bits per spike.
SInfoperspike = SInfo/meanRate;
end