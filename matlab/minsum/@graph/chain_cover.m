function [n_st] = chain_cover(G)
	nV = get(G,'n');
	nE = get(G,'m');
	E = G.E;
	%n_st = network_flow();
		n_st = sdpvar(nE,1);
		
	in = G.in;
	for i=1:length(in)
		in{i}(2,:) = i;
	end
	in{1} = [];
	in{end} = [];
	in = [in{:}];
	
	out = G.out;
	for i=1:length(out)
		out{i}(2,:) = i;
	end
	out{1} = [];
	out{end} = [];	
	out = [out{:}];	
		
	i1 = in(2,:); j1 = in(1,:);
	i2 = out(2,:); j2 = out(1,:);
	A_eq = sparse([i1(:); i2(:)],[j1(:); j2(:)],[ones(1,length(i1)) -ones(1,length(i2))],nV,nE,length(i1)+length(i2));
	f = zeros(nE,1);
	f([G.out{1}])=1;
	b_eq = zeros(nV,1);
	lb = ones(nE,1);
	
	%ops = optimset('LargeScale','on','Display','iter');
	ops = optimset('LargeScale','on','Display','off');
	x = linprog(f,[],[],A_eq,b_eq,lb,inf(nE,1),[],ops);
	n_st = x;
end