function G = sp_grid(G,M,di)

sz = size(M);
n1 = sz(1);
n2 = sz(2);


nV = nnz(M);
G.V = 1:nV;

E = {};
nE = 0;
si = [];

[i j] = find(M);
M = sparse(i,j,G.V,n1,n2,nV);
[i j v] = find(M);

for k=1:size(di,2)
	ii2 = msub2ind(sz,i+di(1,k),j+di(2,k));
	m = zeros(length(ii2),1);
	m(~isnan(ii2)) = M(ii2(~isnan(ii2)));

	E{1,k} = [v(m>0) m(m>0)]';
	si = [si k*ones(1,size(E{1,k},2))];
end


G.structure = di;
G.structure_ind = si;
G.E = cell2mat(E);

G = edge_index(G);

end