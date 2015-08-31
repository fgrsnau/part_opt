function z = mbin2dec(x)
if isempty(x)
	z = 0;
else
	if ischar(x)
		%z = bin2dec(x);
		pows = length(x)-find(x=='1')+1;
		z = sum(bitset(uint32(0),pows));
	else
		pows = length(x)-find(x)+1;
		z = sum(bitset(uint32(0),pows));
		%z = bin2dec(num2str(x));
	end
end

end