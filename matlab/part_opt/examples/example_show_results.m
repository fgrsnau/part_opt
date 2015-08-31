% show segmented image:
figure(1); clf; imagesc(reshape(x,sz),[0 K]); axis equal; axis off;
J = flipud(gray(K+1)); J(1,:) = [1 0 0]; 
colormap(J);

%% TRW-S history:
figure(2); clf;
plot(stats.tt,stats.tLB,'-b','LineWidth',2); hold on; 
plot(stats.tt,stats.tE,'-r','LineWidth',2);
legend({'LB','E'}); xlabel('time, s.'); ylabel('Energy');
set(gcf,'Name','TRW-S history');

%% Partial optimality:
figure(4); clf;
% mask where exactly one label remains alive
mask = sum(X,1)==1;
imagesc(reshape(x.*mask, sz),[0 K]); axis equal; axis off;
colormap(J);
cb = colorbar;
set(gcf,'Name','Proved Optimal part');
set(cb,'YLim',[0 K]);
set(cb,'YTick',[0:K]);
tt = {'unknown'};
for k=1:K
	tt = [tt num2str(k)];
end
set(cb,'YTickLabel',tt);

figure(5); clf;
% number of remaining alive labels
Xn = sum(X,1);
imagesc(reshape(Xn,sz),[0 K]); axis equal; axis off;
colormap(J); cb = colorbar;
set(gcf,'Name','Problem Reminder');
set(cb,'YLim',[0 K]);
set(cb,'YTick',[0:K]);