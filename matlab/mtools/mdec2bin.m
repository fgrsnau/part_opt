function x = mdec2bin(z,d)
	%z = uint32(z);
	x = logical(bitget(z,d:-1:1));
	% x = dec2bin(z,d)=='1';
end