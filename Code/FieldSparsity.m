function sparsityIndex = FieldSparsity(t, o)
%Computes the selectivity of a tuning curve t given an occupancy o as:
%1 - sparsity index. The sparsity index is defined as by Jung, Wiener and 
%McNaughton (JNS, 1994). Closer to one means sparser, so more selective

sparsityIndex = 1 - nansum(t(:).*o(:)/nansum(o(:))).^2/nansum(t(:).^2.*o(:)/nansum(o(:)));

end