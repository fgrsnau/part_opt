function E1 = window(E,W)
%
% E1 = windows(E,W)
% Input:
%		W \subset V -- window
%
%
%
W = row(W);
nV = E.G.get_nV();
nE = E.G.get_nE();
G1 = graph();
G1.V = W;
%G1.W = W;
%G1.V = 1:numel(W);
%mask = sparse(1,W,true,1,numel(E.G.V),numel(W));
ind = sparse(1,W,1:numel(W),1,nV,numel(W));
inner = ind>0;
ee = find(inner(E.G.E(1,:)) & inner(E.G.E(2,:)));
G1.E = full(ind(E.G.E(:,ee)));
G1 = G1.edge_index();
% calculate boundary
B = logical(sparse(1,nV));
st = 1:nE;
outbnd  = inner(E.G.E(1,:)) & ~inner(E.G.E(2,:));
inbnd = ~inner(E.G.E(1,:)) & inner(E.G.E(2,:));
B(E.G.E(1,outbnd)) = true;
B(E.G.E(2,inbnd)) = true;

G1.B = false(size(W));
G1.B(ind(B)) = true;
E1 = energy(G1,E.K);
%E1.B = find(B);
E1.f0 = 0;
E1.f1 = E.f1(:,W);
f2 = get_f2_full(E);
E1.f2 = f2(:,:,ee);
%
%
%
end