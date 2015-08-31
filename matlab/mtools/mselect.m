function B = mselect(A,i1,varargin)
%
%		B = mselect(A,i1,i2..,in, val)
%

switch length(varargin)
	case 0
        B = A(i1);
	case 1
		i2 = varargin{1};
		if(length(i2)==1)
			i2 = repmat(i2,size(i1));
		end
		if(length(i1)==1)
			i1 = repmat(i1,size(i2));
		end
		B = A(sub2ind(size(A),i1,i2));
	case 2
		i2 = varargin{1};
		i3 = varargin{2};
		sz1 = size(i1);
		if numel(i1)==1
			sz1 = 1;
		end
		sz2 = size(i2);
		if numel(i2)==1;
			sz2 = 1;
		end
		sz3 = size(i3);
		if numel(i3)==1;
			sz3 = 1;
		end
		SZ = max(max(sz1,sz2),sz3);
		if numel(i1)==1
			i1 = repmat(i1,SZ);
		end
		if numel(i2)==1
			i2 = repmat(i2,SZ);
		end
		if numel(i3)==1
			i3 = repmat(i3,SZ);
		end
		sz = [length(i1) length(i2) length(i3)];
		if min(sz)==0
			B = zeros([1 0]);
			return;
		end
		if(length(i3)==1)
			i3 = repmat(i3,max(sz));
		end
		B = A(sub2ind(size(A),i1,i2,i3));
  otherwise
        error('selectiong for 4 and more subindeces is not implemented');
end

end