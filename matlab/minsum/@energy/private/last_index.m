function x = last_index(x,i)

switch(ndims(x))

	case 1
		x = x(i);
	case 2
		x = x(:,i);
	case 3
		x = x(:,:,i);
	case 4
		x = x(:,:,:,i);
	otherwise
		error('This large dimension is not implemented');
end

end