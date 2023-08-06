function SInfoperspike = SpatialInfo(t, o)
%Computes the spatial information in tuning curve t, considering an
%occupancy o. This is returned in bits per spike.

%Mean rate of the cell
meanRate = nansum(t(:) .* o(:) / nansum(o(:)));

%Spatial information in bits per seconds.
SInfo = nansum(o(:) / nansum(o(:)) .* t(:) .* log2(t(:) / meanRate));

%Converting in bits per spike.
SInfoperspike = SInfo/meanRate;
end