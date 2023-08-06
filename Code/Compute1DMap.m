function map = Compute1DMap(Xd, Z, nbins) %[mapModel, mapOcc] = ComputeMaps(Xd, Z, params)
% map = Compute1DMap(Xd, Z, nbins)
% Compute 1D map efficiently. 
%
% INPUTS:
% - Xd: the binned independent variable
% - Z: the dependent variable.
% - nbins: scalar, number of bins in Xd over which to accumulate values of 
%   Z.
%
% OUTPUT:
% - map: a 1 x nbins array of Z values summed into vbins of Xd.
%
% USAGE:
% map = Compute1DMap(Xd, Z, nbins);
%
% J. Fournier 07 / 2023

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