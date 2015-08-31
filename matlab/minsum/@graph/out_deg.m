function d = out_deg(G)

d = full(sum(G.out>0,2));

end