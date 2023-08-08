function map = Compute2DMap(Xd, Yd, Z, nXbins, nYbins)
% Compute2DMap - Compute 2D map efficiently by accumulating Z into values of Xd and Yd.
%
%
%   Compute2DMap efficiently computes a two-dimensional map by accumulating
%   values of the dependent variable Z into the bins of the binned independent
%   variables Xd and Yd. It utilizes sparse matrix operations for faster
%   computation.
%
%   map = Compute2DMap(Xd, Yd, Z, nXbins, nYbins)
%
%   INPUTS:
%   Xd:         The binned independent variable along the X-axis.
%   Yd:         The binned independent variable along the Y-axis.
%   Z:          The dependent variable to be accumulated into Xd and Yd bins.
%   nXbins:     Scalar, number of bins in Xd for accumulating Z values.
%   nYbins:     Scalar, number of bins in Yd for accumulating Z values.
%
%   OUTPUT:
%   map:        An nYbins x nXbins array of Z values summed into bins of Xd and Yd.
%
%   USAGE:
%   map = Compute2DMap(Xd, Yd, Z, nXbins, nYbins);
%
%
%   SEE ALSO:
%   Compute1DMap, GaussianSmooth1D, GaussianSmooth, MapsAnalyses1D,
%   MapsAnalyses2D
%
%   Written by J. Fournier in August 2023 for the iBio Summer school.


%%

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