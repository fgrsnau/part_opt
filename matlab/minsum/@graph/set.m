function G = set(G,varargin)
%
%		set(G,propName, propVal)
%
%		propName is one of:
%			'V' -- set of verticies
%			'E'	-- set of edges
%

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
   case 'V'
      G.V = val;
   case 'E'
      G.E = val;
   otherwise
      error('Invalid property')
   end
end