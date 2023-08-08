function DirectionalityIndex = FieldDirectionality(t1, t2)
% FieldDirectionality Compute the directionality index between two tuning curves.
%
% DirectionalityIndex = FieldDirectionality(t1, t2) computes the directionality
% index between two tuning curves t1 and t2. The directionality index measures
% the discrepancy between the tuning curves.
%
% INPUTS:
% - t1: First tuning curve.
% - t2: Second tuning curve.
%
% OUTPUT:
% - DirectionalityIndex: Directionality index between the two tuning curves.
%
% USAGE:
% DirectionalityIndex = FieldDirectionality(t1, t2);
%
% Written by J.Fournier 08/2023 for the iBio Summer school


DirectionalityIndex = abs(sum(t1(:) - t2(:)) ./ sum(t1(:) + t2(:)));

end