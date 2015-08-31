%% Example with a stereo model

commandwindow;
run ../../mtools/mpath.m;
maddpath('../../minsum/');
maddpath('../');

% Load images
I1 = im2double(imread('tsukuba2.bmp'));
I2 = im2double(imread('tsukuba1.bmp'));

figure(1); clf;
imshow([I1 I2]);

K = 15; %number of labels
N1 = size(I1,1);
N2 = size(I1,2)-K; % guarantee N2+K is insied the image
block = 2; %window size for matching in a smaller resolution
n1 = floor(N1/block); 
n2 = floor(N2/block);

sz = [n1 n2];
G = graph(sz); %grid graph of size n1 x n2
nV = G.get_nV(); % number of verticies in the graph
nE = G.get_nE(); % number of edges

%energy function
E = energy(G,K);

tic
% local costs calculation
f1 = zeros(K,nV);
for k = 1:K
	dI = sum((I1(1:N1,1:N2,:)-I2(1:N1,[1:N2]+k,:)).^2,3); % pixelwise SSD
	dI(:,end+1)=0;
	dI(end+1,:)=0;
	II = cumsum(cumsum(dI,1),2); 
	IIsp = II([0:n1]*block+1,[0:n2]*block+1);
	IIq = diff(diff(IIsp(),1,1),1,2);
	f1(k,:) = floor(IIq(:)*1000); % block qualities for disparity k
end

% pairwise cost
for k1=1:K
    for k2=1:K
		d = abs(k1-k2);
		% here we pick a weight resulting in a non-zero integrality
        g_tt(k1,k2) = min(d,3)*10;
    end;
end;
f2 = repmat(g_tt,[1 1 nE]);

%set costs to energy
E.f1 = f1;
E = E.set_f2_full(f2);
fprintf('data preparation %3.2f s.\n',toc);

%% call partial optimality method (doc part_opt_TRWS)
tic
% set the sensetivity parameter to allow method tolerate numerical / input
% cost inaccuracy
ops = struct('sensetivity',-0.01);
[x X stats] = part_opt_TRWS(E,[],[],[],ops);
toc;
fprintf('minimization %3.2f s. (%3.2fs.)\n',toc,stats.time);

%% display results
fprintf('TRW-S bound: %d\n',stats.LB);
fprintf('TRW-S solution: %d\n',stats.E);
example_show_results;
