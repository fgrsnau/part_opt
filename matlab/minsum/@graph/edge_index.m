function G = edge_index(G)

nV = length(G.V);
nE = size(G.E,2);

G.out = sparse(G.E(1,:),[1:nE],[1:nE],nV,nE,nE);
G.in = sparse(G.E(2,:),[1:nE],[1:nE],nV,nE,nE);

G.outT = G.out';
G.inT = G.in';

% G.out = cell(length(G.V),1);
% 	G.in = cell(length(G.V),1);
% 	
% 	for i=1:get(G,'m')
% 		v1 = G.E(1,i);
% 		v2 = G.E(2,i);
% 		G.out{v1} = [G.out{v1} i];
% 		G.in{v2} = [G.in{v2} i];
% 	end
end