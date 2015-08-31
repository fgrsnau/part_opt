function mu_st = get_mu2_index(L)
K = L.K;

nV = get(L.G,'n');
nE = get(L.G,'m');

n = nE*(2*K)+nV+nE;
m = nE*K^2+nV*K;

mu_st = reshape([nV*K+1:m],K,K,[]); %index for pairwise mus

end