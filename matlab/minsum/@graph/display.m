function display(A)
%
% Description:
% ------------
%
%    Displays short info on the object.
%

fprintf(['class <a href="matlab: help ' class(A) ';">' class(A) '</a>{\n']);
%fprintf(['class <a href="matlab: a = ' class(A) '; display(a)' ';">' class(A) '</a>{\n']);
fprintf('\n');
ff = fieldnames(A);
for i=1:length(ff)
	ll(i) = length(ff{i});
end
W = max(ll);
for i=1:length(ff)
	p = repmat(' ',1,W-length(ff{i}));
	fprintf('\t');
	fprintf(p);
	fprintf(ff{i});
	%fprintf(['<a href="matlab: display(' 'A.' ff{i} ');">' ff{i} '</a>']);
	fprintf(': ');
	a = eval(['A.' ff{i}]);
	sz = size(a);
	if isnumeric(a) && prod(sz) == 1
		fprintf(num2str(a));
	else
		fprintf('[');
		fprintf('%gx',sz(1:end-1));
		fprintf('%g ',sz(end));
		if isobject(a)
			fprintf(['<a href="matlab: help ' class(a) ';">' class(a) '</a>']);
			%fprintf(['<a href="matlab: a = ' class(a) '; display(a)' ';">' class(a) '</a>']);
		else
			fprintf(class(a));
		end
		if issparse(a)
			fprintf(' sparse');
		end
		fprintf(']');
	end
	fprintf('\n');
end
fprintf('\n');
mm1 = methods(A);
mm = {};k=1;
for i=1:length(mm1)
	name = mm1{i};
	if strcmp(name,'display') || strcmp(name,'subsasgn') || strcmp(name,'subsref')
		continue;
	end
	mm{k} = name;
	k = k+1;
end

m_disp = mm;
for i=1:length(mm)
	if strcmp(mm{i},'operator_func')
		m_disp{i} = 'operator()';
	end
end

for i=1:length(m_disp)
	ll(i) = length(m_disp{i});
end
W = max(ll);
N = length(mm);
columns = ceil((N-1)/4+1);

% precompute column widths
w = zeros(columns);
for j = 1:4 %% 4 methods in a column, j - row
	for i=1:ceil((N-1)/4+1) %% i - column
		r = (i-1)*4+j; % method number
		if(r<=N)
			w(i) = max(w(i),length(mm{r}));
		end
	end
end

for j = 1:4 %% 4 methods in a column, j - row
	fprintf('\t');
	for i=1:columns %% i - column
		r = (i-1)*4+j; % method number
		if(r<=N)
			s = m_disp{r};
			p = repmat(' ',1,w(i)-length(s)+1);
			fprintf(['<a href="matlab: help ' class(A) '/' mm{r} ';">' s '</a>' p]);
			fprintf([' <a href="matlab: edit ' class(A) '/' mm{r} ';">[.m]</a>  ']);
		end
	end
	fprintf('\n');
end
fprintf('}\n');
end