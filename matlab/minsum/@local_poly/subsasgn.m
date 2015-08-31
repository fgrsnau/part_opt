function A = subsasgn(A,S,B)
%
%

if(length(S)==1)
 	if( isequal(S(1).type,'.'))
 		eval(['A.' S(1).subs ' = B;']);
	else
 		error('() and {} operators are not defined');
	end
else
 	if( isequal(S(1).type,'.'))
 		eval(['A.' S(1).subs ' = subsasgn(A.' S(1).subs ',S(2:end),B);']);
	else
 		error('() and {} operators are not defined');
	end
 end

end