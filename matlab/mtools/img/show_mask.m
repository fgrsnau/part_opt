function show_mask(mask,c,varargin)
%
%	show_mask(mask,c) -- show the boundary of the mask
%


if(~exist('c','var') || isempty(c))
c = 'r';
end

style = '-';

if ischar(c)
	style = [style c];
end

if(exist('mask','var') && ~isempty(mask))

contour(mask,1,style,varargin{:});
%maskbd = bwboundaries(mask);
%length(maskbd)
%{
for i=1:length(maskbd)
	bd = maskbd{i};
	plot(bd(:,2),bd(:,1),c,'LineWidth',1);
end
%}

end