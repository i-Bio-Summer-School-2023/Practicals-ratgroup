function [amp, maxphs] = fitSin(t)
% fitSin Fits a sinusoid to a tuning curve and computes amplitude and phase.
%
% [amp, maxphs] = fitSin(t) fits a sinusoid to a tuning curve t and computes
% the amplitude (amp) and phase of the sinusoid.
%
% INPUT:
% - t: 1D vector, typically a tuning curve.
%
% OUTPUTS:
% - amp: Amplitude of the best-fitting sinusoid.
% - maxphs: Phase of the peak of the best-fitting sinusoid in degrees.
%
% USAGE:
% [amp, maxphs] = fitSin(t);
%
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"

t = t(:)';  % Ensure t is a row vector
nphsbins = numel(t);
Phi = (0.5:(nphsbins-0.5))/nphsbins*2*pi;

% Numerator and denominator for phase calculation
Phi_num = sum((t-mean(t)).*sin(Phi));
Phi_den = sum((t-mean(t)).*cos(Phi));

%Phase of the best fitting sinusoid in radians.
phs0 = atan2(Phi_den,Phi_num);

%Amplitude of the best fitting sinusoid.
amp = sum((t-mean(t)).*sin(Phi + phs0))/sum(sin(Phi).^2);

%Phase of the peak of the best fitting sinusoid in degrees.
maxphs = mod(90 - 360*phs0/(2*pi),360);
end