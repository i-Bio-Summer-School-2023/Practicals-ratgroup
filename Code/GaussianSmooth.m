function output = GaussianSmooth(input, smthNbins)
%Smooth a nD array (input) with a gaussian of s.d. = smthNbins (in bins
%unit). SmthNbins should have the same number of elements as the number of
%dimension in input.

%Building the gaussian function that'll be used to smooth the data
Ndim = ndims(input);
npts = 5 * ones(1, Ndim);
npts(smthNbins == 0) = 0;
Smoother = ones(round(2 * npts .* max(smthNbins,1) + 1)); 
for k = 1:Ndim
    x = (-(npts(k) * max(1,smthNbins(k))):(npts(k) * max(1,smthNbins(k))));
    if smthNbins(k) > 0
        Smoother_1D = exp(-x.^2 / smthNbins(k)^2 / 2);
    else
        Smoother_1D = double(x == 0);
    end
    Smoother_1D = Smoother_1D(:);
    vperm = 1:Ndim;
    vperm([1 k]) = vperm([k 1]);
    Smoother_1D = permute(Smoother_1D,vperm);
    sz = size(Smoother);
    sz(k) = 1;
    Smoother = Smoother .* repmat(Smoother_1D, sz);
end
Smoother = Smoother / sum(Smoother(:));

%Detecting NaN values to exclude them from the convolution
valididx = ~isnan(input);

%Replacing NaN values with 0 so they don't count when convolving
input(isnan(input)) = 0;

%Convolving the input vector with the gaussian
output = convn(input,Smoother,'same');

%Counting the actual number of valid points smoothed (i.e. excluding NaN
%values.
flat = convn(double(valididx),Smoother,'same');

%Normalizing the convolved vector by the number of valid points
output = output ./ flat;

%Replacing back NaN values in the output vector.
output(~valididx) = NaN;
end