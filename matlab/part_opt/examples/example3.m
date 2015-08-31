%% Example with a model from an OpenGM file

commandwindow;
run ../../mtools/mpath.m;
maddpath('../../minsum/');
maddpath('../');

pth = './'; fname = 'snail.h5';

% We can pass the filename to the solver
tic
ops = struct([]);
[x X stats] = part_opt_TRWS([pth fname],[],[],ops);
toc;
fprintf('minimization %3.2f s. (%3.2fs.)\n',toc,stats.time);
% try to guess the image size of the model for displaying
sz = size(auto_reshape(x));
% display results
fprintf('TRW-S bound: %d\n',stats.LB);
fprintf('TRW-S solution: %d\n',stats.E);
example_show_results;

return
%% Execute this cell to continue

% % We can also load the modet into matlab
[ee f1 f2] = opengm_read_mex([pth '/' fname],'gm');
% And construct @energy class from it
% create graph from the list of edges in ee
K = size(f1,1);
nV = size(f1,2);
nE = size(f2,3);
G = graph();
G.V = 1:nV;
G.E = double(ee+1);
G = G.edge_index(); % init graph
E = energy(G,K);
E = E.set_f1(f1);
E = E.set_f2_full(f2);

% And evaluate
E.cost(x)
% Or run another method
[X elim P time] = invoke_po_method('Kovtun',E);
[ans x] = max(X,[],1); x = x.*(sum(X)==1); % partial labeling
fprintf('Kovtun time: %3.2fs.\n',time);
example_show_results;
