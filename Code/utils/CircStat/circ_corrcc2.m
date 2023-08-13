function [rho,pval] = circ_corrcc2(alpha1, alpha2, w)
%
% [rho pval ts] = circ_corrcc(alpha1, alpha2)
%   Circular correlation coefficient for two circular random variables.
%
%   Input:
%     alpha1	sample of angles in radians
%     alpha2	sample of angles in radians
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%   Topics in circular statistics, S.R. Jammalamadaka et al., p. 176
%
% PHB 6/7/2008
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

% compute mean directions
n = sum(w(:));
alpha1_bar = circ_mean(alpha1(:), sum(w,2));
alpha2_bar = circ_mean(alpha2(:), (sum(w,1))');

alpha1 = repmat(alpha1(:),[1 size(w,2)]);
alpha2 = repmat(alpha2(:)',[size(w,1) 1]);
% compute correlation coeffcient from p. 176
num = sum(w(:).*(sin(alpha1(:) - alpha1_bar) .* sin(alpha2(:) - alpha2_bar)));
den = sqrt(sum(w(:).*sin(alpha1(:) - alpha1_bar).^2) .* sum(w(:).*sin(alpha2(:) - alpha2_bar).^2));
rho = num / den;	

% compute pvalue
l20 = sum(w(:).*sin(alpha1(:) - alpha1_bar).^2)/sum(w(:));
l02 = sum(w(:).*sin(alpha2(:) - alpha2_bar).^2)/sum(w(:));
l22 = sum(w(:).*(sin(alpha1(:) - alpha1_bar).^2) .*(sin(alpha2(:) - alpha2_bar).^2))/sum(w(:));

ts = sqrt((n * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));