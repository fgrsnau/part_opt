function E = reparam(E,M)
%
%	E = reparam(E,M) -- add reparametrization M
%	M [K x nE x 2] double

E = f2_to_full(E);
G = E.G;
nE = G.get_nE();
nV = G.get_nV();
K = E.K;

f1 = E.f1;
f2 = get_f2_full(E);

f2 = f2 - reshape(kron(M(:,:,2),ones(E.K,1)),K,K,[]);
f2 = f2 - reshape(kron(ones(E.K,1),M(:,:,1)),K,K,[]);

i1 = kron(G.E(2,:)',ones(K,1));
f1 = f1 + accumarray([repmat((1:K)',nE,1) i1], reshape(M(:,:,2),[],1),[K nV]);

i1 = kron(G.E(1,:)',ones(K,1));
f1 = f1 + accumarray([repmat((1:K)',nE,1) i1], reshape(M(:,:,1),[],1),[K nV]);

m1 = min(f1,[],1);
%m1(:) = 0;
f1 = f1-repmat(m1,K,1);

m2 = min(reshape(f2,K*K,[]),[],1);
f2 = f2 - reshape(repmat(m2,K*K,1),K,K,[]);

E.f1 = f1;
E = set_f2_full(E,f2);
E.f0 = E.f0+sum(m1)+sum(m2);

end