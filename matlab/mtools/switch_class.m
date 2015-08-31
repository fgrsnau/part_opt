function switch_class(var,class)
		evalin('caller',[var '=' class '(' var ')']);
end