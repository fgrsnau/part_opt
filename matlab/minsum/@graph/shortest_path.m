function [S val] = shortest_path(G,s,t,wV,wE)
% s,t - verticies
%	wE - edge cost matrix	
%

nV = get_nV(G);

dist = inf(1,nV);
prev = nan(1,nV);
dist(s) = wV(s);

open = logical(sparse(nV,1));
open(s) = true;

while any(open)
	vv = find(open);
	[ans i] = sort(dist(vv),'ascend');
	u = vv(i(1));
	open(u) = false;
	if(u==t)% read path
		val = dist(u);
		S = [u];
		while ~isnan(prev(u))
			u = prev(u);
			S = [u S];
		end
	end
	out = find(G.out(u,:));
	for e = out
		v = G.E(2,e);
		d = dist(u)+wE(e)+wV(v);
		if(d<dist(v))
			dist(v) = d;
			prev(v) = u;
			open(v) = true;
		end
	end
end
end