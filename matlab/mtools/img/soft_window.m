function H = soft_window(sz)
	h1 = 1/2 - 1/2*cos(2*pi*[1:sz(1)]'/sz(1));
	h2 = 1/2 - 1/2*cos(2*pi*[1:sz(2)]/sz(2));
	H = h1*h2;
end
