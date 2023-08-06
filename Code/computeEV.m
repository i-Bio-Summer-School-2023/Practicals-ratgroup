function EV = computeEV(y, ypred)
%

RSS = nansum((y - ypred).^2);
m = nanmean(y);
TSS = nansum((y - m).^2);

EV = 1 - (RSS/TSS);
end