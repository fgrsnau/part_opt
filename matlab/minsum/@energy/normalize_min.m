function E = normalize_min(E)

E = f2_to_full(E);

G = E.G;
nE = G.get_nE();
nV = G.get_nV();
K = E.K;

f1 = E.f1;
f2 = get_f2_full(E);

m1 = min(f1,[],1);
f1 = f1-repmat(m1,K,1);

m2 = min(reshape(f2,K*K,[]),[],1);
f2 = f2 - reshape(repmat(m2,K*K,1),K,K,[]);

E.f1 = f1;
E = set_f2_full(E,f2);
E.f0 = E.f0+sum(m1)+sum(m2);

end