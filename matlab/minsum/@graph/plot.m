function plot(G,sz)

set(gcf,'Color',[1 1 1]);

i1 = 1+mod(G.V-1,sz(1));
i2 = 1+floor((G.V-1)/sz(1));
hh = ishold();
%plot(i2,i1,'ok','MarkerSize',5); hold on;
%plot([reshape(i2(G.E(1,:)),1,[]); reshape(i2(G.E(2,:)),1,[])],[reshape(i1(G.E(1,:)),1,[]); reshape(i1(G.E(2,:)),1,[])],'-k','LineWidth',2);
X1 = reshape(i2(G.E(1,:)),1,[]);
X2 = reshape(i2(G.E(2,:)),1,[]);
Y1 = reshape(i1(G.E(1,:)),1,[]);
Y2 = reshape(i1(G.E(2,:)),1,[]);

q = quiver(X1,Y1,X2-X1,Y2-Y1,0,'-ok');
set(q,'MaxHeadSize',0.15,'ShowArrowHead','off','MarkerFaceColor',[1 1 1],'Marker','o','MarkerSize',5);
set(q,'ShowArrowHead','on');
%get(q,'ShowArrowHead')
%get(q,'MaxHeadSize')
%refreshdata(q);
%,'-k','LineWidth',2);
if hh
	hold on;
else
	hold off;
end
axis off;
axis ij; axis equal;
end