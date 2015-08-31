function neib = get_neib(G,v)
% neib = get_in(G,v)	- get neighboring verticies
% 
	neib = [G.E(1,get_in(G,v)) G.E(2,get_out(G,v))];
end