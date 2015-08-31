%% This example will segment random noisy image with Potts model segmentation 
%% into K labels

commandwindow;
run ../../mtools/mpath.m;
maddpath('../../minsum/');
maddpath('../');

%% Create a grid graph
N = 100; sz = [N N];
G = graph(sz); % grid graph N x N nodes
K = 5; E = energy(G,K); % Energy over graph G with 5 labels

nV = get(G,'n'); % number of vertices in the graph
nE = get(G,'m'); % number of edges in the graph

rand('seed',1);randn('seed',1);
E = set_f1(E,rand(K,nV)); % set univariate potentials randomly

f2_pots = repmat(0.2*(1-eye(K,K)),[1 1 nE]); % create Potts model pairwise potentials
E = set_f2_full(E,f2_pots); % set to energy

%% call partial optimality method (doc part_opt_TRWS)
[x X stats] = part_opt_TRWS(E);
% Can calculate energy of x:
fprintf('Resulting energy: %d\n',cost(E,x));
%% display results
fprintf('TRW-S bound: %d\n',stats.LB);
fprintf('TRW-S solution: %d\n',stats.E);
example_show_results;

return

%% compare to a reference solution
commandwindow;
[x2 LB] = alg_trws(E,0.00001,200);
fprintf('TRW-S bound: %d\n',LB);
fprintf('TRW-S solution: %d\n',cost(E,x2));
