function pth = file_path(f)

	[pathstr, name, ext] = fileparts(f);
	pth =  pathstr;

end