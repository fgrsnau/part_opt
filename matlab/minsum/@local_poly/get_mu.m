function [mu1 mu2] = get_mu(L,x)
%
%		[mu1 mu2] = get_mu(L,x)
%
%		Input:
%			x : [1xnV] -- configuration
%   Output:
%			mu1 : [KxnV] -- node-marginals
%			mu2 : [KxKxnE] -- edge-marginals


nV = get(L.G,'n');
nE = get(L.G,'m');
K = L.K;
EE = get(L.G,'E');

mu1 = sparse(x(:),1:nV,ones(nV,1),K,nV);

i1 = x(EE(1,:));
i2 = x(EE(2,:));
i3 = 1:nE;

mu2 = zeros(K,K,nE);
mu2 = write(mu2,i1,i2,i3,ones(1,nE));

end