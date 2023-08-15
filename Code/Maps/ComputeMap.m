function map = ComputeMap(Xd, Yd, Z, nXbins, nYbins)
% ComputeMap - Compute a map efficiently by accumulating Z into values of Xd and Yd.
%
%   map = ComputeMap(Xd, Yd, Z, nXbins, nYbins) efficiently computes a two-dimensional
%   map by accumulating values of the dependent variable Z into the bins of the binned
%   independent variables Xd and Yd. It utilizes sparse matrix operations for faster
%   computation.
%
% INPUTS:
% - Xd: The binned independent variable along the X-axis. If empty, computes a 1D map along Y.
% - Yd: The binned independent variable along the Y-axis. If empty, computes a 1D map along X.
% - Z: The dependent variable to be accumulated into Xd and Yd bins.
% - nXbins: Scalar, number of bins in Xd for accumulating Z values.
% - nYbins: Scalar, number of bins in Yd for accumulating Z values.
%
% OUTPUT:
% - map: An nYbins x nXbins array of Z values summed into bins of Xd and Yd.
%
% USAGE:
% map = ComputeMap(Xd, Yd, Z, nXbins, nYbins);
%
% SEE ALSO:
% GaussianSmooth, MapsAnalyses
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological recordings"
%%

%If Yd is empty, we will compute a 1D map along X
if isempty(Yd)
    Yd = ones(size(Xd));
    nYbins = 1;
end

%If Xd is empty, we will compute a 1D map along Y
if isempty(Xd)
    Xd = ones(size(Yd));
    nXbins = 1;
end

%Selecting valid indices
valididx = ~isnan(Xd) & ~isnan(Yd) & ~isnan(Z);
Xd = Xd(valididx);
Yd = Yd(valididx);
Z = Z(valididx);

%Summing Z within indices of Xd.
map = sparse(Yd, Xd, Z, nYbins, nXbins);

%Converting into a full matrix again for accessing elements more
%conveniently
map = full(map);

end