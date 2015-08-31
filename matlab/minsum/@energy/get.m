function val = get(E,propName)
%
%		get(G,prop), where prop is one of:
%			'G' -- graph            [graph]
%			'K'	-- number of labels
%			'f1' -- univariate func [K x nV] 
%			'f2' -- pairwise func   [K x K x nE]
%
%		See also GRAPH
%

switch propName
case 'G'
   val = E.G;
case 'K'
   val = E.K;
case 'f1'
   val = E.f1;
case 'f2'
   error('depricated. f2 should be accessed using dedicated methods.');
case 'L'
	val = E.L;
otherwise
   error([propName ,' is not a valid property'])
end