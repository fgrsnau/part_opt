function E = set_f2(E,f2,edges)
%
% E = set_f2(E,f2) sets sparse pairwise potentials to given list of edges
%
%	Input:
%		E - [@energy]
%		f2 - [@user_class] -- pairwise potentials
%													must implement
%															f_cost(f2,x_st)
%															f_matrix(f2)
%															f_min_sum_conv(f2,v)
%															f_tmin_sum_conv(f2,v)
%   edges - list of edges
%
%

	ind = length(E.f2.store{4}.data)+1;
	E.f2.store{4}.data{ind} = f2;
	E.f2.index.type(edges) = 4;
	E.f2.index.i(edges) = ind;

end
