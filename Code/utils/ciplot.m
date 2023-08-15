function plot_handle = ciplot(x,lower,upper,colour,alpha)
% ciplot(x,lower,upper)
% ciplot(x,lower,upper,colour)
% ciplot(x,lower,upper,colour,alpha)
%
% Plots a shaded region on a graph between specified lower and upper confidence
% intervals.
%
% INPUT:
% - x: X data values (optional, default to index values if empty).
% - lower: Lower confidence interval values (vector).
% - upper: Upper confidence interval values (vector).
% - colour: Color of the shaded region (optional, default is [0.5 0.5 0.5]).
% - alpha: Transparency of the shaded region (optional, default is 0.5).
%
% OUTPUT:
% - plot_handle: Handle to the filled area plot.
%
% USAGE:
% ciplot(x, lower, upper, colour, alpha);
%
% Note:
% - lower and upper vectors must be of the same length.
% - Multiple shaded plots can be overlaid without problems by using 'fill' 
%   instead of 'area' and making them transparent.
%
% originally created by Raymond Reynolds 24/11/06
% modified by:
% Pham Thai Binh 12/06/2017
% J. Fournier 08/2023

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end
if nargin<5
    alpha=0.5;
end
if nargin<4
    colour= [0.5 0.5 0.5];
end
if isempty(x)
    x=1:length(lower);
end
% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
    x=x'; 
end
if find(size(lower)==(max(size(lower))))<2
    lower=lower'; 
end
if find(size(upper)==(max(size(upper))))<2
    upper=upper'; 
end
if nargout == 0
    fill([x fliplr(x)],[upper fliplr(lower)], colour, 'EdgeColor', 'none', 'FaceAlpha',alpha);
else
    plot_handle = fill([x fliplr(x)],[upper fliplr(lower)], colour, 'EdgeColor', 'none', 'FaceAlpha',alpha);
end
