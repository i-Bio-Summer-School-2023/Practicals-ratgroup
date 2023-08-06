function map = Compute2DMap(Xd, Yd, Z, nXbins, nYbins)
%Compute 2D map efficiently. Xd and Yd are the binned independent variables
%and Z is the dependent variable. nXbins and nYbins correspond to the
%number of discretizing bins in Xd and Yd respectively.

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