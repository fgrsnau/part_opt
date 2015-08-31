function out = get_out(G,W)
% G.get_out(v) -- edges going out of set W
%out = G.out{v};
if numel(W)==1
    %out = find(G.out(W,:));
	%out = nonzeros(G.out(W,:));
	out = find(G.outT(:,W))';
else
    %[ans1 ans2 out] = find(G.out(W,:));
	[ans1 ans2 out] = find(G.outT(:,W));
	out = col(out);
end
end