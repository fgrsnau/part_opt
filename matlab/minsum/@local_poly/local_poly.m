function L = local_poly(varargin)
% 	local_poly(G,K) -- create a local polytope.
% 	local_poly() -- create empty local polytope.
%

switch (nargin)
	case 0
		L = [];
		L.G = [];
		L.K = 0;
	case 2
		L = [];
		L.G = varargin{1};
		L.K = varargin{2};
	otherwise
		error('wrong arguments');		
end

L = class(L,'local_poly');

end

