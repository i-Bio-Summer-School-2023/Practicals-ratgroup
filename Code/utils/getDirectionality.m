function DirectionalityIndex = getDirectionality(t1, t2)
% getDirectionality Compute the directionality index between two tuning curves.
%
% DirectionalityIndex = getDirectionality(t1, t2) computes the directionality
% index between two tuning curves t1 and t2. The directionality index measures
% the discrepancy between the tuning curves as the absolute difference
% divided by the sum
%
% INPUTS:
% - t1: First tuning curve.
% - t2: Second tuning curve.
%
% OUTPUT:
% - DirectionalityIndex: Directionality index between the two tuning curves.
%
% USAGE:
% DirectionalityIndex = getDirectionality(t1, t2);
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%

DirectionalityIndex = abs(sum(t1(:) - t2(:)) ./ sum(t1(:) + t2(:)));

end