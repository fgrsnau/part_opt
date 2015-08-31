function X = auto_reshape(x)
n = length(x);
s = sqrt(n);
for i = floor(s/2): 2*floor(s) % go over reasonable form factors
	if mod(n,i)==0 % divisible
		X = reshape(x,i,n/i);
		dx = diff(X,1,1);
		%dy = diff(X,1,2);
		dxy = diff(dx,1,2);
		C(i) = sum(sum(abs(dxy)));
		%fprintf('%i x %i score: %f\n',i,n/i,C(i));
		%cfigure(1);
		%clf; plot(C,'.r'); drawnow;
	else
		C(i) = inf;
	end
end
C(C==0) = inf;
[v i] = min(C);
X = reshape(x,i,n/i);
end