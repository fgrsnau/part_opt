function f2 = get_f2_full(E)
%
if strcmp(E.type,'full')
	f2 = E.f2;
	return;
end
%
nV = get(E.G,'n');
nE = get(E.G,'m');
K = E.K;

f2 = zeros(K,K,nE);
tt = unique(E.f2.index.type);
if(isempty(tt))
	return
end
for t = tt
	if t==1
		ii = E.f2.index.type==t;
		%index = E.f2.index.i(ii);
        r = E.f2.store{t}.data;
        if(size(r,3)==1)
            r = repmat(r,[1 1 length(ii)]);
        end
		f2(:,:,ii) = r;
	else
		ii = E.f2.index.type==t;
		index = E.f2.index.i(ii);
		if(iscell(E.f2.store{t}.data))
			for i = 1:length(ii)
				if(ii(i))
					f2(:,:,i) = f2(:,:,i) + f2_matrix(E,i);%E.f2.store{t}.f_matrix(E.f2.store{t}.data{index(i)});
				end
			end
		elseif t==1
			f2(:,:,ii) = f2(:,:,ii) + E.f2.store{t}.f_matrix(E.f2.store{t}.data(:,:,index));
		else
			f2(:,:,ii) = f2(:,:,ii) + E.f2.store{t}.f_matrix(last_index(E.f2.store{t}.data,index));
		end
	end
end

end