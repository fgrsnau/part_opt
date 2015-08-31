function E = apply_mapping(E,P)
%
%	E = apply_mapping(E,M) -- apply permutation (mapping) P
%	P [K x nV] double

P = double(P);

E = f2_to_full(E);
G = E.G;
nE = G.get_nE();
nV = G.get_nV();
K = E.K;

if ndims(P)==2
	%f1 = E.f1;
	%f2 = get_f2_full(E);
	[i s] = mmeshgrid([K,nV]);
	i1 = mselect(P,i,s);
	E.f1 = E.f1-mselect(E.f1,i1,s);
	st = repmat(reshape(1:nE,1,1,nE),[K,K,1]);
	i = repmat(reshape(1:K,K,1,1),[1,K,nE]);
	j = repmat(reshape(1:K,1,K,1),[K,1,nE]);
	s = mselect(G.E,1,st);
	t = mselect(G.E,2,st);
	i1 = mselect(P,i,s);
	j1 = mselect(P,j,t);
	E.f2 = E.f2-mselect(E.f2,i1,j1,st);
	E.f0 = 0;
	%
else
	for s=1:nV
		g1 = zeros(K,1);
		for i1=1:K
			for i=1:K
				g1(i1) = g1(i1) + E.f1(i,s)*P(i,i1,s);
			end
		end
		E.f1(:,s) = E.f1(:,s) - g1;
	end
	for e=1:nE
		s = G.E(1,e);
		t = G.E(2,e);
		g2 = zeros(K,K);
		for i1=1:K
			for j1=1:K
				for i=1:K
					for j=1:K
						g2(i1,j1) = g2(i1,j1) + E.f2(i,j,e)*P(i,i1,s)*P(j,j1,t);
					end
				end
			end
		end
		E.f2(:,:,e) = E.f2(:,:,e) - g2;
	end
	E.f0 = 0;
end
%
end