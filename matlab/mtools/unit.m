function a = unit(x)

if islogical(x)
	a = 1;
else
	a = ones(1,1,class(x));
end

end