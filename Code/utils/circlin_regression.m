function [rho ,pval, s ,phi0] = circlin_regression(X, Phi)
% Calculates the linear-circular coefficient (rho) between a linear variable X
% and a circular variable Phi (in radians). Also returns the corresponding
% p-value, the slope (s), and the intercept (phi0) of the best fitting line.
% See Kempter et al., 2012.
%
% INPUT:
% - X: Linear variable (can contain NaN values).
% - Phi: Circular variable in radians (can contain NaN values).
%
% OUTPUT:
% - rho: Linear-circular coefficient.
% - pval: p-value of the significance test for rho.
% - s: Slope of the best fitting line.
% - phi0: Intercept of the best fitting line.
%
% USAGE:
% [rho, pval, s, phi0] = circlin_regression(X, Phi);
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"

% Remove NaN values from X and Phi
valididx = ~isnan(X) & ~isnan(Phi);
X = X(valididx);
Phi = Phi(valididx);

n = numel(X);
nresol = 3601; % Resolution for angular range
ncycles = 50; % Number of cycles in angular range
A = linspace(-ncycles,ncycles,nresol); % Range of angular frequencies

% Calculate the resultant length R for different angular frequencies (A)
R = NaN(nresol,1);
for i = 1:numel(A)
    a = A(i);
    R(i) = sqrt((1/n * sum(cos(Phi - 2*pi*a*X)))^2 + (1/n * sum(sin(Phi - 2*pi*a*X)))^2);
end

% Find the angular frequency (slope s) that maximizes the resultant length 
% R
[~,imax] = max(R);
s = A(imax);

% Calculate phi0 (intercept of the best fitting line)
phi0_num = sum(sin(Phi - 2*pi*s*X));
phi0_den = sum(cos(Phi - 2*pi*s*X));
phi0 = atan2(phi0_num,phi0_den);

% Compute the circular mean of X transformed to angular data
theta = mod(2*pi*abs(s)*X, 2*pi);
theta_num = sum(sin(theta));
theta_den = sum(cos(theta));
theta_bar = atan2(theta_num,theta_den);

% Compute the circular mean of Phi
Phi_num = sum(sin(Phi));
Phi_den = sum(cos(Phi));
Phi_bar = atan2(Phi_num,Phi_den);

% Calculate the circular correlation coefficient (rho)
rho_num = sum(sin(Phi - Phi_bar).*sin(theta - theta_bar));
rho_den = sqrt(sum(sin(Phi - Phi_bar).^2)*sum(sin(theta - theta_bar).^2));
rho = rho_num/rho_den;

% Calculate test statistics and p-value for significance
lambda22 = 1/n * sum((sin(Phi - Phi_bar).^2).*(sin(theta - theta_bar).^2));
lambda02 = 1/n * sum((sin(Phi - Phi_bar).^0).*(sin(theta - theta_bar).^2));
lambda20 = 1/n * sum((sin(Phi - Phi_bar).^2).*(sin(theta - theta_bar).^0));
z = rho * sqrt(n * (lambda02 * lambda20) / lambda22);
pval = 1 - erf(abs(z)/sqrt(2));
end