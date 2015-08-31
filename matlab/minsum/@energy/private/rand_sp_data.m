function D = rand_sp_data(sz)

D.support = logical(sprand(sz(1),sz(2),0.3));
[D.i1 D.i2] = find(D.support);
D.support = sparse(D.i1,D.i2,[1:length(D.i1)],sz(1),sz(2));
D.cost = rand(1,length(D.i1));
D.max_val = max(D.cost)+rand(1);

end