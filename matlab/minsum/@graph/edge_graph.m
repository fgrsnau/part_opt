function B = edge_graph(A)
nE = get_nE(A);
B = graph();
B.V = 1:nE;
%how many 
nBE = sum(in_deg(A).*out_deg(A));
B.E = zeros(2,nBE);
k=1;
for e1 =1:nE
	ee = find(A.out(A.E(2,e1),:));
	for e2 = ee
		B.E(:,k) = [e1,e2]';
		k = k+1;
	end
end

B = edge_index(B);

end