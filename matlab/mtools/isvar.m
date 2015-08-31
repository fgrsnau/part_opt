function a = isvar(Name)
	a = evalin('caller',['exist(''' Name ''',''var'')']);
end