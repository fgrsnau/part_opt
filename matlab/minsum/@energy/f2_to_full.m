function E = f2_to_full(E)
%
%	E = f2_to_full(E) -- converts all special matricies to full matricies
%

if(~strcmp(E.type,'full')) % && any(E.f2.index.type>1))
	f2 = get_f2_full(E);
	%E = set_f2_full(E,f2);
	E.f2 = f2;
	E.type = 'full';
end

end