function [rho ,pval, s ,phi0] = circlin_regression(X, Phi)
%Calculates the linear-circular coefficient (rho) between linear variable X
%and circular variable Phi (in radians). Also returns the corresponding 
%p-value, the slope and the intercept of the best fitting line. See Kempter
%et al, 2012.

valididx = ~isnan(X) & ~isnan(Phi);
X = X(valididx);
Phi = Phi(valididx);

n = numel(X);
nresol = 3601;
ncycles = 50;
R = NaN(nresol,1);
A = linspace(-ncycles,ncycles,nresol);
for i = 1:numel(A)
    a = A(i);
    R(i) = sqrt((1/n * sum(cos(Phi - 2*pi*a*X)))^2 + (1/n * sum(sin(Phi - 2*pi*a*X)))^2);
end

[~,imax] = max(R);
s = A(imax);

phi0_num = sum(sin(Phi - 2*pi*s*X));
phi0_den = sum(cos(Phi - 2*pi*s*X));
phi0 = atan2(phi0_num,phi0_den);

theta = mod(2*pi*abs(s)*X, 2*pi);
theta_num = sum(sin(theta));
theta_den = sum(cos(theta));
theta_bar = atan2(theta_num,theta_den);

Phi_num = sum(sin(Phi));
Phi_den = sum(cos(Phi));
Phi_bar = atan2(Phi_num,Phi_den);

rho_num = sum(sin(Phi - Phi_bar).*sin(theta - theta_bar));
rho_den = sqrt(sum(sin(Phi - Phi_bar).^2)*sum(sin(theta - theta_bar).^2));
rho = rho_num/rho_den;

lambda22 = 1/n * sum((sin(Phi - Phi_bar).^2).*(sin(theta - theta_bar).^2));
lambda02 = 1/n * sum((sin(Phi - Phi_bar).^0).*(sin(theta - theta_bar).^2));
lambda20 = 1/n * sum((sin(Phi - Phi_bar).^2).*(sin(theta - theta_bar).^0));
z = rho * sqrt(n * (lambda02 * lambda20) / lambda22);

pval = 1 - erf(abs(z)/sqrt(2));
end