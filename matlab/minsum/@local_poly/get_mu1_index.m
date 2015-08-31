function mu_s = get_mu1_index(L)
K = L.K;

nV = get(L.G,'n');
nE = get(L.G,'m');

n = nE*(2*K)+nV+nE;
m = nE*K^2+nV*K;

mu_s = reshape([1:nV*K],K,[]);  %index for univariate mus

end