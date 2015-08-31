function [XX elim P time] = invoke_po_method(method,E,X0,y,eps0,UB)

K = E.K;
G = E.G;
nV = E.get_nV();
nE = E.get_nE();
if(~isvar('X0'))
	X0 = ones([K,nV],'int32');
end

if strcmpi('Iterative_mex_TRWS',method)
	tic;
	%y = [];
	ops = struct([]);
	%ops = struct('it_batch',1);
	[y1 P M time hist burn] = part_opt_TRWS_mex(int32(G.E-1),double(E.f1),double(E.f2),int32(X0),[],int32(y-1),ops);
	y1 = double(y1)+1;
	%assert(all(y1==y));
	t1 = toc;
	% do something non-parallel
	%{
	fprintf('-------saving---------\n');
	%cfigure(1); clf(); drawnow;
	f = fopen('1.tmp','w');
	for i=1:1000
		fprintf(f,[datestr(now) '\n']);
	end
	fclose(f);
	%}
	%save('test','M');
	%clear part_opt_TRWS_mex;
	P = P+1;
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
end

if strcmpi('Swoboda0_TRWS',method)
	tic;
	[P M time] = part_opt_TRWS_mex(int32(G.E-1),int32(E.f1),int32(E.f2),int32(X0),[],int32(y-1),double(1));
	t1 = toc;
	clear part_opt_TRWS_mex;
	P = P+1;
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
end

if strcmpi('PBP_TRWS',method)
	tic;
	[X P M time] = PBP_mex(int32(G.E-1),int32(E.f1),int32(E.f2),int32(X0),[],int32(y-1));
	t1 = toc;
	clear PBP_mex;
	P = P+1;
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
end

if strcmpi('L1',method)
	[X elim P time] = po_L1(E,X0,y,eps0);
	XX = P == repmat([1:K]',1,nV);
	%elim = nnz(XX<X0);
end
if strcmpi('L1_UB',method)
	if ~exist('UB','var') || isempty(UB)
		UB = E(y);
	end
	fprintf('UB=%f\n',UB);
	[X elim P time] = po_L1_ub(E,X0,y,eps0,UB);
	XX = P == repmat([1:K]',1,nV);
	%elim = nnz(XX<X0);
end

if strcmpi('Iterative_CPLEX',method)
	[XX elim P time] = po_swoboda_ip2(E,X0,y,eps0);
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
end
if strcmpi('Swoboda0_CPLEX',method)
	[XX elim P time] = po_swoboda_ip0(E,X0,y);
	%po_swoboda_ip2(E,X0,y,eps0);
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
end

if strcmpi('Swoboda_CPLEX',method)
	tic;
	[XX elim P] = po_swoboda_ip(E,X0,y);
	time = toc;
end

if strcmpi('Kovtun',method)
	%tic;
	[XX elim P time] = po_kovtun_plain(E,X0,[]);
	%time = toc;
end

if strcmpi('MQPBO',method)
	%tic;
	[LB XX P time] = po_MQPBO(E,X0,0);
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
	%time = toc;
end

if strcmpi('MQPBO-P',method)
	%tic;
	[LB XX P time] = po_MQPBO(E,X0,1);
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
	%time = toc;
end

if strcmpi('DEE',method)
	tic;
	[LB XX P] = po_dee(E,X0,1);
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
	time = toc;
end

if strcmpi('DEE2',method)
	tic;
	[LB XX P] = po_dee2(E,X0,1);
	XX = P == repmat([1:K]',1,nV);
	elim = nnz(XX<X0);
	time = toc;
end

end