function [amp, maxphs] = fitSin(t)
%Computes the amplitude and phase of the sinusoid that best fits tuning
%curve t.

t = t(:)';
nphsbins = numel(t);
Phi = (0.5:(nphsbins-0.5))/nphsbins*2*pi;

Phi_num = sum((t-mean(t)).*sin(Phi));
Phi_den = sum((t-mean(t)).*cos(Phi));

%Phase of the best fitting sinusoid in radians.
phs0 = atan2(Phi_den,Phi_num);

%Amplitude of the best fitting sinusoid.
amp = sum((t-mean(t)).*sin(Phi + phs0))/sum(sin(Phi).^2);

%Phase of the peak of the best fitting sinusoid in degrees.
maxphs = mod(90 - 360*phs0/(2*pi),360);
end