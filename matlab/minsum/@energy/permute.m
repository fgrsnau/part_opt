function E = permute(E,P)
%
%	E = reparam(E,M) -- apply permutation (mapping) P
%	P [K x nV] double

E = f2_to_full(E);
G = E.G;
nE = G.get_nE();
nV = G.get_nV();
K = E.K;

%f1 = E.f1;
%f2 = get_f2_full(E);
[i s] = mmeshgrid(K,nV);
i1 = mselect(P,i,s);
E.f1 = mselect(E.f1,i1,s);
st = repmat(reshape(1:nE,1,1,nE),K,K,1);
i = repmat(reshape(1:K,K,1,1),1,K,nE);
j = repmat(reshape(1:K,1,K,1),K,1,nE);
s = G.E(1,st);
t = G.E(2,st);
i1 = mselect(P,i,s);
j1 = mselect(P,j,t);
E.f2 = mselect(E.f2,i1,j1,st);

%
%
end