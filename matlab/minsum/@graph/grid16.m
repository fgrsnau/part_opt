function G = grid8(G,sz)

G = grid_reg(G,sz,[1 0;0 1;1 1;-1 1; 2 1; 1 2; -1 2; -2 1]');

end