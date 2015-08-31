function E = add_global_minl(E,support,L1,D1,L2,D2)
%
%  E = add_global_minl(E,L1,D1,L2,D2) 
%  adds the cost min(L1*x+d1,L2*x+D2) to the binary energy function by
%  convering this cost to quadratic submodular using the thechnique by
%  Kohli-Ladicky-Torr
%
%  Input:	E -energy
%			support [nV nf] logical sparse
%			L1, L2 - linear terms over support set support
%			D1, D2 - constants
%
%

%L1 = L1(:)';
%L2 = L2(:)';

if(E.K~=2)
	error('Energy must be in 2 labels.');
end

nf = size(support,2);

if (size(L1,2)~=nf)
	error('L1 must be nV x %i sparse',nf);
end

if (size(L2,2)~=nf)
	error('L2 must be nV x %i sparse',nf);
end

%% check coefficients are of definite signe and arrange needed sign
nnV = [];
for i=1:nf
	dL = L1(support(:,i),i)-L2(support(:,i),i);
	if(any(dL<0) && any(dL>0))
		error('Can only represent when L1-L2 is of definite sign.');
	end
	if(any(dL>0)) %% swap L1(:,i) and L2(:,i)
		R = L2(:,i);
		L2(:,i) = L1(:,i);
		L1(:,i) = R;
		D = D1(i);
		D1(i) = D2(i);
		D2(i) = D;
	end
	nnV(i) = nnz(support(:,i));
end

f2 = get_f2_full(E);

dL = L1-L2;
dD = D1-D2;

%% now represent min(dL*x+dD,0) = min_y y(dL*x+dD) = min_y y*dL*x + y*dD
nV = get_nV(E);
y = nV+1:nV+nf;
E.G.V = [E.G.V y]; % auxiliari verticies

% new edges between y and x from support set
newedges = zeros(2,nnz(support));
k = 1;
for i=1:nf
	newedges(1,k:k+nnV(i)-1) = y(i);
	newedges(2,k:k+nnV(i)-1) = E.G.V(support(:,i));
	k = k+nnV(i);
end
E.G.E = [E.G.E newedges];
E.G = E.G.edge_index(); % update incidence matrix

%% add L2 to linear terms on x
E.f1(2,:) = E.f1(2,:)+sum(L2,2)';

%% add D2 to constant term
E.f0 = E.f0+sum(D2);

%% set linear term dD*y
f1y = zeros(2,nf);
f1y(2,:) = dD;
E.f1 = [E.f1 f1y];

%% set quadratic term y*dL*x
f2y = zeros(2,2,nnz(support));
f2y(2,2,:) = dL(support);
f2 = cat(3,f2,f2y);
E = set_f2_full(E,f2);
%% done
end