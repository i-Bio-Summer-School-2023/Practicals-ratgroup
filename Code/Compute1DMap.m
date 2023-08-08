function map = Compute1DMap(Xd, Z, nbins) %[mapModel, mapOcc] = ComputeMaps(Xd, Z, params)
% Compute1DMap - Compute 1D map efficiently.
%
% Compute1DMap takes a binned independent variable Xd and a corresponding
% dependent variable Z, then it efficiently computes a one-dimensional map
% by accumulating values of Z into the specified number of bins (nbins) of Xd.
% The function utilizes sparse matrix operations for efficient computation.
%
% INPUTS:
%   Xd: The binned independent variable.
%   Z: The dependent variable.
%   nbins: Scalar, number of bins in Xd over which to accumulate values of Z.
%
% OUTPUT:
%   map: A 1 x nbins array of Z values summed into bins of Xd.
%
% USAGE:
%   map = Compute1DMap(Xd, Z, nbins);
%
% SEE ALSO:
%   Compute2DMap, GaussianSmooth1D, GaussianSmooth, MapsAnalyses1D,
%   MapsAnalyses2D
%
% Written by J. Fournier in August 2023 for the iBio Summer school.


%%

%Selecting valid indices
valididx = ~isnan(Xd) & ~isnan(Z);
Xd = Xd(valididx);
Z = Z(valididx);

%Summing Z within indices of Xd.
map = sparse(ones(size(Xd)), Xd, Z, 1, nbins);

%Converting into a full matrix again for accessing elements more
%conveniently
map = full(map);

end