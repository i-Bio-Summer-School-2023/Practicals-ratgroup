function selectivityIndex = FieldSelectivity(t)
%Computes the selectivity of the tuning curve t as the maximal amplitude
%normalized by the mean.

selectivityIndex = (max(t(:)) - min(t(:))) / nanmean(t(:));

end