function c = f2_cost(E,st,labeling)
	t = E.f2.index.type(st);
	i = E.f2.index.i(st);
	if(iscell(E.f2.store{t}.data))
	  if t==4
	    c = E.f2.store{t}.data{i}.f_cost(labeling);
	  else
			c = E.f2.store{t}.f_cost(E.f2.store{t}.data{i},labeling);
		end
	elseif t==1
		c = E.f2.store{t}.f_cost(E.f2.store{t}.data(:,:,i),labeling);
	else
		c = E.f2.store{t}.f_cost(last_index(E.f2.store{t}.data,i),labeling);
	end
end