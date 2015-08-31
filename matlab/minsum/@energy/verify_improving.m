function [v x1 mu1 mu2 time] = verify_improving(E,P,UB)
%
%	E = reparam(E,M) -- apply permutation (mapping) P
%	P [K x nV] double
%   P [K x K x nV] double -- relaxed mapping

if exist('UB','var') && ~isempty(UB)
	f_UB = get_theta(E);
	f_UB = [E.f0-UB f_UB];
else
	f_UB = [];
end
%f_UB = f_UB
E = apply_mapping(E,P);
[x1 v mu1 mu2 phi duals time] = alg_lp(E,f_UB);
%
%
%
end