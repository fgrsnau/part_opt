function G = grid_reg(G,sz,di)
%
%		create regular grid graph with givven connectivity pattern
%		sz -- size of grid
%   di -- [2 x n] double -- list of displacements from each vertex to be connected by an edge
%

n1 = sz(1);
n2 = sz(2);
nV = n1*n2;
G.V = 1:nV;

E = {};
nE = 0;
si = [];

for i=1:size(di,2)
	[i2 i1] = meshgrid(1:n2,1:n1);
	ii1 = msub2ind(sz,i1,i2);
	[i2 i1] = meshgrid(1+di(2,i):n2+di(2,i),1+di(1,i):n1+di(1,i));
	ii2 = msub2ind(sz,i1,i2);
	%m = i1>=0 & i1<sz(1) & i2>=0 & i2<sz(2);
	m = find(~isnan(ii2));
	u = ii1(m);
	v = ii2(m);
	E{1,i} = [u(:) v(:)]';
	si = [si i*ones(1,length(m))];
end


G.structure = di;
G.structure_ind = si;
G.E = cell2mat(E);

G = edge_index(G);

end