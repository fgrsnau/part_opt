function [G ee] = add_edges(G,GE)
%
%  Add edges GE to graph G avoiding duplicate edges, edge index ee will return new indices of GE  
%

[c ia ib] = intersect(G.E',GE','rows');
last_edge = size(G.E,2);
maskb = true(1,size(GE,2));
maskb(ib) = false;

G.E = [G.E GE(:,maskb)];
ee = nan(size(maskb));
ee(maskb) = last_edge+1:last_edge+nnz(maskb);
ee(ib) = ia;

G = edge_index(G);

end