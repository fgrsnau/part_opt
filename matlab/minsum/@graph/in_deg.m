function d = in_deg(G)

d = full(sum(G.in>0,2));

end