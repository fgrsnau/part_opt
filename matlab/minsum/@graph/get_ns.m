function ns = get_ns(G,s)
%
%		ns = get_ns(graph) -- number of edges adjucent to s
%

ns = length(get_in(G,s))+length(get_out(G,s));

end