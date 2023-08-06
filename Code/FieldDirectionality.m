function DirectionalityIndex = FieldDirectionality(t1, t2)
%Computes the discrepancy between to tuning curves t1 and t2.

DirectionalityIndex = abs(sum(t1(:) - t2(:)) ./ sum(t1(:) + t2(:)));

end