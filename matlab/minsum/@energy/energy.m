function E = energy(varargin)
%
% 	energy(G,K) -- create an energy minimization problem instance eover graph G with K labels.
% 	energy() -- empty energy object.
%
%   See also	ENERGY/GET, ENERGY/COST, ENERGY/GET_POLY,
%							ENERGY/GET_THETA, ENERGY/RAND_INSTANCE, ENERGY/SET_F2_FULL
%
%		For full list of methods see: methods(energy)
%

E=[];
E.G = graph();
E.K = 0;
E.f0 = 0;
E.f1 = [];
E.f2 = [];
E.type = 'full';

switch(nargin)
	case 0
		E.L = local_poly();
	case 2
		G = varargin{1};
		K = varargin{2};
		E.G = G;
		E.K = K;
		E.L = local_poly(G,K);
	otherwise
		error('wrong arguments');
end

E = class(E,'energy');
%E = construct_f2(E);
end