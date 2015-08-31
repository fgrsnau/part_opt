function in = get_in(G,v)
% G.get_in(v) -- edges going into v
	%in = G.in{v};
	%in = find(G.in(v,:));
	in = find(G.inT(:,v))';
end