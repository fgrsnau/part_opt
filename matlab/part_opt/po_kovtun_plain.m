function [XX elim P time] = po_kovtun_plain(E,XX,y)
%
% X [KnV x K] int32 - alive nodes
%
K = E.K;
nV = E.get_nV();
f1 = E.f1;
f2 = E.get_f2_full();
ee = E.G.E;
fprintf('____________Kovtun-plain_________________\n');
%
XX0 = XX;
P = int32(repmat([1:K]',1,nV));
time = 0;
for k=1:K
	y1 = ones(1,nV)*k;
	[XX1 P1 st] = part_opt_IK_mex(int32(ee-1),single(f1),single(f2),int32(XX0),int32(y1-1));
	time = time + st.time;
	elim = nnz(XX1<XX0);
	P1 = P1+1;
	P = compose(P1,P);
	fprintf('label%i/%i  eliminated %i\n',k,K,elim);
end

XX = P == repmat([1:K]',1,nV);
elim = nnz(XX<XX0 & f1<2^29);
maxelim = sum(f1(:)<2^29)-nV;
fprintf('eliminated %i (%3.2f%%)\n',elim,elim/maxelim*100);

end