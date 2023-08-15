function output = filterellip(input, low, high, samplingRate)
% filterellip applies elliptic bandpass filtering to the input signal.
%
% output = filterellip(input, low, high, samplingRate) filters the input signal using
% an elliptic bandpass filter with specified low and high cutoff frequencies.
%
% INPUTS:
% - input: Input signal to be filtered.
% - low: Lower cutoff frequency of the bandpass filter (in Hz).
% - high: Upper cutoff frequency of the bandpass filter (in Hz).
% - samplingRate: Sampling rate of the input signal (in Hz).
%
% OUTPUT:
% - output: Filtered output signal.
%
% USAGE:
% output = filterellip(input, low, high, samplingRate);
%
% See also: ellip, filtfilt.
%%
[b,a] = ellip(2, 0.1, 40, [low high] / (samplingRate/2));
input(isnan(input)) = 0;
output = filtfilt(b, a, input);

end