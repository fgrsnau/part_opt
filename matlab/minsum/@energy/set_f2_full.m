function E = set_f2_full(E,f2,edges)
%
% E = set_f2_full(E,f2,edges=[]) sets all (or edges) pairwise potentials in the form of one 3d array
%
%	Input:
%		E - [@energy]
%		f2 - [K x K x nE] -- full 3d array of pairwise potentials
%
%

if strcmp(E.type,'full')
	if exist('edges','var')
		E.f2(:,:,edges) = f2;
	else
		E.f2 = f2;
	end
	return
end

if ~exist('edges','var') || isempty(edges)
	nE = get(E.G,'m');
	K = E.K;
	E.f2.index.i = [1:nE];
	E.f2.index.type = ones(1,nE);
	E.f2.store{1}.data = f2;
else
	if(size(f2,3)==1)
		f2 = repmat(f2,[1,1,length(edges)]);
	end
	E.f2.store{1}.data(:,:,edges) = f2;
	E.f2.index.type(edges) = 1;
	E.f2.index.i(edges) = edges;
end

end