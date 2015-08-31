function E = construct_f2(E)
nV = get(E.G,'n');
nE = get(E.G,'m');
K = E.K;

E.f2 = [];
E.f2.index = [];
E.f2.index.i = ones(1,nE);
E.f2.index.type = ones(1,nE);

E.f2.store = {};
E.f2.store{1} = [];
E.f2.store{1}.type = 'full';
E.f2.store{1}.f_cost = @f2_full_cost;
E.f2.store{1}.f_matrix = @f2_full_matrix;
E.f2.store{1}.f_min_sum_conv = @f2_full_min_sum_conv;
E.f2.store{1}.f_tmin_sum_conv = @f2_full_tmin_sum_conv;

E.f2.store{2} = [];
E.f2.store{2}.type = 'potts';
E.f2.store{2}.f_cost = @f2_potts_cost;
E.f2.store{2}.f_matrix = @f2_potts_matrix;
E.f2.store{2}.f_min_sum_conv = @f2_potts_min_sum_conv;
E.f2.store{2}.f_tmin_sum_conv = @f2_potts_min_sum_conv;

E.f2.store{1}.data = zeros(K,K,1);
E.f2.store{2}.data = zeros(K,1);


E.f2.store{3} = [];
E.f2.store{3}.type = 'sparse';
E.f2.store{3}.f_cost = @f2_sparse_cost;
E.f2.store{3}.f_matrix = @f2_sparse_matrix;
E.f2.store{3}.f_min_sum_conv = @f2_sparse_min_sum_conv;
E.f2.store{3}.f_tmin_sum_conv = @f2_sparse_tmin_sum_conv;

E.f2.store{3}.data = {};
E.f2.store{3}.data{1} = {};
E.f2.store{3}.data{1}.support  = sparse([]);
E.f2.store{3}.data{1}.i1 = [];
E.f2.store{3}.data{1}.i2 = [];
E.f2.store{3}.data{1}.cost = [];
E.f2.store{3}.data{1}.max_val = 0;


E.f2.store{4} = [];
E.f2.store{4}.type = 'user_class';
E.f2.store{4}.data = {};
end

function C = f2_full_matrix(D)
	C = D;
end
function C = f2_full_cost(D,x_st)
	C = D(x_st(1),x_st(2));
end
function y = f2_full_min_sum_conv(D,x)
y = min(D+repmat(x',size(D,1),1),[],2);
end
function y = f2_full_tmin_sum_conv(D,x)
y = min(D+repmat(x,1,size(D,2)),[],1)';
end

function c = f2_potts_cost(dd,x_st)
if x_st(1)==x_st(2)
	c = dd(x_st(1));
else
	c = 0;
end
end

function f2 = f2_potts_matrix(dd)
f2 = reshape(full(spdiags(dd,size(dd,1)*[0:size(dd,2)-1],size(dd,1),size(dd,1)*size(dd,2))),size(dd,1),size(dd,1),[]);
end

function y = f2_potts_min_sum_conv(d,x)
y = min(min(x),x+d);
end


function C = f2_sparse_cost(D,x_st)
  i = D.support(x_st(1),x_st(2));
  if i>0
     C = D.cost(i);
  else
     C = D.max_val;
  end
end

function f2 = f2_sparse_matrix(D)
  f2 = accumarray([D.i1(:) D.i2(:)],D.cost,size(D.support),[],D.max_val);
end

function y = f2_sparse_min_sum_conv(D,x)
v = D.cost(:)+x(D.i2);
y = accumarray(D.i1(:),v,[size(D.support,1) 1],@min,inf);
y = min(y,min(x)+D.max_val);
end

function y = f2_sparse_tmin_sum_conv(D,x)
v = D.cost(:)+x(D.i1);
y = accumarray(D.i2(:),v,[size(D.support,2) 1],@min,inf);
y = min(y,min(x)+D.max_val);
end