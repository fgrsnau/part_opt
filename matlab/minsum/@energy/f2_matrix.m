function c = f2_matrix(E,st)
	t = E.f2.index.type(st);
	i = E.f2.index.i(st);
	if(iscell(E.f2.store{t}.data))
	  if t==4
	    c = E.f2.store{t}.data{i}.f_matrix();
	  else
		  c = E.f2.store{t}.f_matrix(E.f2.store{t}.data{i});
		end
	else
		c = E.f2.store{t}.f_matrix(last_index(E.f2.store{t}.data,i));
	end
end