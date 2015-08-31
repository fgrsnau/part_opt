function E = random_instance(E,type,varargin)
nV = get(E.G,'n');
nE = get(E.G,'m');
K = E.K;

if(~exist('type','var'))
	type = [];
end

E = construct_f2(E);
%E.type = type;

if strfind(type,'full')
    E = f2_to_full(E);
end

if strfind(type,'potts')
	E.type = 'potts';
	if(nargin>2)
		sigma = varargin{1};
	else
		sigma = 1;
	end

	if strfind(type,'N')
		E.f1 = randn(K,nV);
		E.f2.index.i = [1:nE];
		E.f2.index.type(:) = 2;
		lambda = sigma*randn(1,nE);
		%E.f0 = sum(lambda);
		E.f2.store{2}.data = kron(ones(K,1),-lambda);
		if strfind(type,'-')
			E.f2.store{2}.data = -abs(E.f2.store{2}.data);
		end
	else
		sign = -1;
		if strfind(type,'+')
			sign = +1;
		end
		E.f1 = rand(K,nV);
		E.f2.index.i = [1:nE];
		E.f2.index.type(:) = 2;
        if strfind(type,'upotts')
            E.f2.store{2}.data = 0.5*repmat(rand(1,nE),K,1)*sign;
        else
            E.f2.store{2}.data = 0.5*rand(K,nE)*sign;
        end
	end

	if strfind(type,'int32')
		E.f1 = floor(100*E.f1);
		E.f2.store{2}.data = floor(100*E.f2.store{2}.data);
	end
	
elseif strfind(type,'sparse')
	E.type = 'sparse';
	E.f1 = rand(K,nV);
	E.f2.index.i = [1:nE];
	E.f2.index.type(:) = 3;
	E.f2.store{3}.data = {};
	for st=1:E.G.get_nE()
		E.f2.store{3}.data{st} = rand_sp_data([K K]);
		if strfind(type,'int32')
			E.f2.store{3}.data{st}.cost = floor(100*E.f2.store{3}.data{st}.cost);
			E.f2.store{3}.data{st}.max_val = floor(100*E.f2.store{3}.data{st}.max_val);
		end
	end
	
else %default% if strfind(type,'full')
	if strcmp(E.type,'full')
		E.type = 'full';
		if strfind(type,'int32')
			E.f1 = int32(floor(100*rand(K,nV)));
			E.f2 = int32(floor(100*rand(K,K,nE)));
		else
			E.f1 = floor(100*rand(K,nV));
			E.f2 = floor(100*rand(K,K,nE));
		end
	else
		E.f1 = rand(K,nV);
		E.f2.index.i = [1:nE];
		E.f2.index.type(:) = 1;
		E.f2.store{1}.data = rand(K,K,nE);
		if strfind(type,'int32')
			E.f1 = floor(100*E.f1);
			E.f2.store{1}.data = floor(100*E.f2.store{1}.data);
		end
	end
end

end