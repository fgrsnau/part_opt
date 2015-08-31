function G = graph(varargin)
% Create a graph.
% 	graph() -- empty graph
% 	graph([n1 n2]) -- grid graph n1 x n2, 4-connected
%
%
G = [];
G.V = [];
G.E = zeros(2,0);
G.out = [];
G.in = [];
G.outT = [];
G.inT = [];
G.B = [];
G.structure = [];
G.structure_ind = [];
switch nargin
	case 0
	case 1
		%% grid graph
		sz = varargin{1};
		n1 = sz(1);
		n2 = sz(2);
		nV = n1*n2;
		G.V = 1:nV;
		nE = (n1-1)*n2+(n2-1)*n1;
		%% create edge -> vertex index
		[i2 i1] = meshgrid(1:n2-1,1:n1);
		ii1 = sub2ind(sz,i1,i2);
		[i2 i1] = meshgrid(2:n2,1:n1);
		ii2 = sub2ind(sz,i1,i2);
		Eh = [ii1(:)'; ii2(:)'];
		[i2 i1] = meshgrid(1:n2,1:n1-1);
		ii1 = sub2ind(sz,i1,i2);
		[i2 i1] = meshgrid(1:n2,2:n1);
		ii2 = sub2ind(sz,i1,i2);
		Ev = [ii1(:)'; ii2(:)'];
		G.E = [Eh Ev];
	otherwise
		error('Unrecognized arguments');
end
G = class(G,'graph');
G = edge_index(G);
end
