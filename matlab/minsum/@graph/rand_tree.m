function G = rand_tree(G)
%
%	Random spanning tree of G. This is bad O(VE) algorithm
%

nV = get_nV(G);
G.E = [G.E flipdim(G.E,1)]; % add reverse edges

E = zeros(2,nV-1);
T = zeros(1,nV);
v = randi(nV);
T(v) = 1; % mark node is in the tree
for i = 1:nV-1
	ee = find(T(G.E(1,:)) & ~T(G.E(2,:)));
	e = ee(randi(length(ee)));
	E(:,i) = G.E(:,e);
	T(G.E(:,e)) = 1;
end

G.E = flipdim(flipdim(E,1),2);
G = edge_index(G);

end