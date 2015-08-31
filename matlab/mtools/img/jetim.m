function I = jetim(A,c,palette)
if ~exist('c','var') || isempty(c)
	a = min(A(:));
	b = max(A(:));
	%m1 = min(a,-b);
	%m2 = max(b,-a);
	m1 = a;
	m2 = b;
else
	m1 = c(1);
	m2 = c(2);
	A = max(A,m1);
	A = min(A,m2);
end
if ~exist('palette','var') || isempty(palette)
	palette = jet;
end

Ai = floor((A-m1)/(m2-m1)*(size(palette,1)-1))+1;
mask = ~isnan(Ai);
I = nan(numel(A),3);
I(mask,:) = palette(Ai(mask),:);
%I(~mask,:) = nan;
I = reshape(I,[size(A,1) size(A,2) 3]);

end