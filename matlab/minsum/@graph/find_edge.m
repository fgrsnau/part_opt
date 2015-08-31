function e = find_edge(G,st)

s = st(1);
t = st(2);

o = find(G.out(s,:));
i = find(G.E(2,o)==t);
e = o(i);

end