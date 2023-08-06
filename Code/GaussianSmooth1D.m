function output = GaussianSmooth1D(input, smthNbins)
% output = GaussianSmooth1D(input, smthNbins)
%
% Smooth a 1D vector with a gaussian.
%
% INPUT:
% - input : 1D array to be smoothed.
% - smthNbins: s.d. of the Gaussian.
% 
% OUTPUT:
% - output: smoothed 1D array of the same size as input.
%
% USAGE:
% output = GaussianSmooth1D(input, smthNbins);
%
% J Fournier 07/2023

%Saving sise of input to reshape it at the end
sz = size(input);
input = input(:);

%Building the gaussian function that'll be used to smooth the data
npts = 5;
x = (-(npts * max(1,smthNbins)):(npts * max(1,smthNbins)));

if smthNbins > 0
    Smoother = exp(-x.^2/smthNbins^2/2);%Gaussian
else
    Smoother = double(x == 0);%dirac
end
Smoother = Smoother(:) / sum(Smoother(:));

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

%Reshaping the output vector as originally.
output = reshape(output, sz);

end