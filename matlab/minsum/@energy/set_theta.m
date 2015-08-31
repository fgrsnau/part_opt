function E = set_theta(E,theta)
%
%		E = set_theta(E,theta) -- sets all potentials
%
%		Input:
%			E [@energy]
%			theta [K*nV + K*K*nE] -- concatenated vector of potentials
%
	l = E.K*E.G.get('n');
	E.f1 = reshape(theta(1:l),E.K,[]);
	f2 = reshape(theta(l+1:end),[E.K, E.K, E.G.get('m')]);
	E = set_f2_full(E,f2);
end