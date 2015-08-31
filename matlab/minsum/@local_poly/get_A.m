function [A_eq b_eq] = get_A(L)
%
%		[A_eq b_eq] = get_A(L) returns equality constraint matrix of local polytope
%

K = L.K;

nV = get(L.G,'n');
nE = get(L.G,'m');

n = nE*(2*K)+nV+nE;
m = nE*K^2+nV*K;

mu_s = get_mu1_index(L);  %index for univariate mus
mu_st = get_mu2_index(L); %index for pairwise mus
%mu_s = reshape([nE*K^2+1:m],K,[]);  %index for univariate mus
%mu_st = reshape([1:nE*K^2],K,K,[]); %index for pairwise mus
mu_ts = mu_st;
for i=1:nE
	mu_ts(:,:,i) = mu_st(:,:,i)';
end

i1 = kron([1:nV],ones(1,K));

A_norm = sparse(i1(:),mu_s(:),ones(nV*K,1),nV,m,nV*K);
b_norm = ones(nV,1);

i1 = kron([1:K*nE],ones(1,K));

i2 = [1:K*nE]; E = get(L.G,'E'); j2 = mu_s(:,E(2,:));
A_marg1 = sparse([i1 i2],[mu_st(:); j2(:)],[ones(1,length(i1)) -ones(1,length(i2))] ,K*nE,m,K*(K+1)*nE);
b_marg1 = zeros(K*nE,1);

i1 = kron([1:K*nE],ones(1,K));
i2 = [1:K*nE]; E = get(L.G,'E'); j2 = mu_s(:,E(1,:));
A_marg2 = sparse([i1 i2],[mu_ts(:); j2(:)],[ones(1,length(i1)) -ones(1,length(i2))] ,K*nE,m,K*(K+1)*nE);
b_marg2 = zeros(K*nE,1);

A_eq = [A_norm; A_marg1; A_marg2];
b_eq = [b_norm; b_marg1; b_marg2];

end