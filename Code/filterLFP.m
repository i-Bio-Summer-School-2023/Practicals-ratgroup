function [filteredLFP] = filterLFP(input, low, high, sampRate)
% Filters input using a elliptic filter in a band
% from 'low' to 'high' Hz. sampRate is the sampling frequency of the input
% signal.

[b,a] = ellip(2, 0.1, 40, [low high]/(sampRate/2));
input(isnan(input)) = 0;
filteredLFP = filtfilt(b, a, input);