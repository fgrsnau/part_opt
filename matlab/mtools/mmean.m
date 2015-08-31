function A = mmean(A, dims)
W = ~isnan(A); % weights
A(~W) = 0;
dims = sort(dims);
for i=length(dims):-1:1
	A = sum(A,dims(i));
	W = sum(W,dims(i));
end
A = A./W;
end