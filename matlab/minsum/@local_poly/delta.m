function mu = delta(L,x)
%
%		mu = delta(L,x)
%
%		Input:
%			x : [1xnV] -- configuration
%   Output:
%			mu : -- linear embeddings of x
%
%

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

mu = [full(mu1(:)); full(mu2(:)); 1];

end