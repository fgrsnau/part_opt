function val = get(G,propName)
%
%		get(G,prop), where prop is one of:
%			'V' -- set of verticies
%			'E'	-- set of edges
%			'n' -- number of vert.
%			'm' -- number of edges.
%
switch propName
case 'E'
   val = G.E;
case 'V'
   val = G.V;
case 'n'
   val = length(G.V);
case 'm'
   val = size(G.E,2);
otherwise
   error([propName ,'Is not a valid property'])
end