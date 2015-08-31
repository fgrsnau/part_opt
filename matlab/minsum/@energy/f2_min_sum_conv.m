function y = f2_min_sum_conv(E,st,x)
	if strcmp(E.type,'full')
		D = double(E.f2(:,:,st));
		y = min(D+repmat(x',size(D,1),1),[],2);
		return;
	end
	t = E.f2.index.type(st);
	i = E.f2.index.i(st);
	if(iscell(E.f2.store{t}.data))
		if(t==4)
			y = f_min_sum_conv(E.f2.store{t}.data{i},x);
		else
			y = E.f2.store{t}.f_min_sum_conv(E.f2.store{t}.data{i},x);
		end
	else
 		y = E.f2.store{t}.f_min_sum_conv(last_index(E.f2.store{t}.data,i),x);
	end
end