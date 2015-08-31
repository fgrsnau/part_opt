function varargout = subsref(A,S)
% Implements operator. to access member variables and methods and
% operator() to call operator_func.m method if exists
%
% Handles arbitrary chain of operators, for example
% f.translate([2 2]).smooth(5).fplot()
% is the same as fplot(smooth(translate(f,[2 2]),5))
%
% ! Does not work from within class method for that class
%
% (C) Alexander Shekhovtsov
%
if isequal(S(1).type,'.') % first operator is '.'
	if isfield(struct(A),S(1).subs) %% member object
		if length(S)>1
			%% perform rest of operators on the member object
			[varargout{1:nargout}] = subsref(eval(['A.' S(1).subs]),S(2:end));
		else
			%% just return member object
			[varargout{1:nargout}] = eval(['A.' S(1).subs]);
		end
	else % not member object
		if length(S)==1   % . and no more operators
			if ismethod(A, S(1).subs) %method but () operator was not specified
				eval(['help ' class(A) '/' S(1).subs]);
				error(['"' S(1).subs '": function call is missing']);
			else % neither member nor method
				error(['"' S(1).subs '" is neither a member nor a method of class ' class(A)]);
			end
		else  % . and more operators and not member object => think it is method
			if length(S)>2
				% call method(A,arguments) and the rest of operators
				[varargout{1:nargout}] = subsref(eval([S(1).subs '(A,S(2).subs{:})']),S(3:end));
			else
				% call method(A,arguments)
				[varargout{1:nargout}] = eval([S(1).subs '(A,S(2).subs{:})']);
			end
		end
	end
elseif isequal(S(1).type,'()') % first operator is ()
	if length(S)>1
		% call operator_func(arguments) and the rest of operators
		[varargout{1:nargout}] = subsref(operator_func(A,S(1).subs{:}),S(2:end));
	else
		% call operator_func(arguments)
		[varargout{1:nargout}] = operator_func(A,S(1).subs{:});
	end
else
	error(['operator {} is not defined for class ' class(A)]);
end
end