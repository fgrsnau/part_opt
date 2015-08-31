function E = set_f2_sparse(E,f2,edges)
%
% E = set_f2_sparse(E,f2,edges) sets sparse pairwise potentials to given list of edges
%
%	Input:
%		E - [@energy]
%	   f2 - [K x K] -- sparse 2D array of pairwise potentials
%   edges - list of edges
%
%

	if (islogical(f2))
		D = {};
		D.support = f2;
		[D.i1 D.i2] = find(D.support);
		D.support = sparse(D.i1,D.i2,[1:length(D.i1)],size(f2,1),size(f2,2));
		D.cost = zeros(1,length(D.i1));
		D.max_val = 1e5;%inf;
	else
		D = {};
		D.support = f2~=0;
		[D.i1 D.i2 D.cost] = find(f2);
		D.support = sparse(D.i1,D.i2,[1:length(D.i1)],size(f2,1),size(f2,2));
		D.max_val = 1e5;%inf;
	end

	ind = length(E.f2.store{3}.data)+1;
	E.f2.store{3}.data{ind} = D;
	E.f2.index.type(edges) = 3;
	E.f2.index.i(edges) = ind;
end
