function E = add_f2_full(E,f2,ee)
%
% E = set_f2_full(E,f2,ee) sets pairwise potentials on edges ee to 3D array f2
%
%	Input:
%		E - [@energy]
%		f2 - [K x K x nee] -- full 3d array of pairwise potentials
%		ee - [nee] -- edge indicies
%


	nE = get(E.G,'m');
	K = E.K;


% edges for which there is already full data -> add to existing

  mask = E.f2.index.type(ee) == 1;
  data_index = E.f2.index.i(ee(mask));
  E.f2.store{1}.data(:,:,data_index) = E.f2.store{1}.data(:,:,data_index)+f2(:,:,mask);


% all remaining edges -> add new data
  f2 = f2(:,:,~mask);
  ee = ee(~mask);
	
	nee = length(ee);
	
	lastfull = size(E.f2.store{1}.data,3);
	
	E.f2.index.i(ee) = [lastfull+1:lastfull+nee];
	E.f2.index.type(ee) = ones(1,nee);
	E.f2.store{1}.data = cat(3,E.f2.store{1}.data,f2);

end