function corr = CalculatePearson(t1,t2,d)
%corr = CalculatePearson(t1,t2,d) returns the Pearson's correlation between
%Nd arrays t1 and t2, integrating over dimensions specified in d. 
%If d is not provided, CalculatePearson operates along the first 
%non-singleton dimension.

sz = size(t1);
if nargin < 3
    d = find(sz > 1, 1, 'first');
end

if ~isequal(sz, size(t2))
    error('arrays should have the same dimension')
end

t1 = t1 - mean(t1,d,'omitnan');
t2 = t2 - mean(t2,d,'omitnan');
corr = sum(t1.*t2, d, 'omitnan')./sqrt(sum(t1.^2, d, 'omitnan').*sum(t2.^2, d, 'omitnan'));

end