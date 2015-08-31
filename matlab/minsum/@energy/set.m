function E = set(E,varargin)
%
%		set(E,propName, propVal)
%
%		propName is one of:
%			'G' -- graph
%			'K'	-- number of labels
%			'f1' -- univariate func
%			'f2' -- pairwise func
%

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
	prop = propertyArgIn{1};
	val = propertyArgIn{2};
	propertyArgIn = propertyArgIn(3:end);
	switch prop
		case 'G'
			E.G = val;
		case 'K'
			E.K = val;
		case 'f1'
			E.f1 = val;
		case 'f2'
			error('depricated. f2 should be accessed using dedicated methods.');
		otherwise
			error([prop ,'Is not a valid property'])
	end
end