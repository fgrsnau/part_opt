function [x X stats] = part_opt_TRWS(varargin)
% [x X stats] = part_opt_TRWS(E,X0,M,y,ops)
% [x X stats] = part_opt_TRWS(file_name,M,y,ops)
% Input 1:
%	0) E -- @energy class
%	3) X0 -- [K x nV] int32 -- unary mask of alive labels (reserved, currently no effect)
%	4) M -- [K x nE] double -- starting backward messages (default [])
%	5) y -- [1 x nV] int32 -- test labeling (default [])
%	6) ops struct double -- options, see available fields from an output
% Input 2:
%	file_name -- opengm model
%	M, y, ops -- as above
%
% Output:
%	x - integer labeling found by TRW-S in the initialization phase
%	X - [K x nV] int32 -- mask of alive labels
% stats - struct various output statistics
%	P - [K x nV] int32 -- imroving mapping
%	M - [K x nE] double -- current backward messages
%	LB - bound from TRW-S
%    E - best energy found in TRW-S
%	time - total runnig time
%
%	tt - time stamps over iterations
%   tLB - history of the lower bound over iterations
%   tE - history of energy over iterations
%   tdM - history of message change over iterations
%
% 	burn [K x nV x nC] - double - progress of prunning:
% 	layer 0: PO iteration number when prunned,
% 	layer 1: margin when prunned
% 	layer 2: prunning condition: { INIT = 0, WTA = 1, CUT = 2, PXCUT = 3}
% 	layer 3: sensetivity  - sensetivity level when pruned (reserved)
% 	layer 4: time when pruned
% 	layer 5: TRW-S iteration number when prunned
% 	*/

if(isa(varargin{1},'energy'))
	E = varargin{1};
	if(nargin>1)
		X0 = varargin{2};
	end
	if(nargin>2)
		M = varargin{3};
	end
	if(nargin>3)
		y = varargin{4};
	end
	if(nargin>4)
		ops = varargin{5};
	end
	
	nV = E.get_nV();
	nE = E.get_nE();
	K = E.K;
	
	if(~isvar('X0') || isempty(X0))
		X0 = ones([1 nE],'int32');
	end
	
	if(~isvar('M'))
		M = [];
	end
	
	if(~isvar('y'))
		y = [];
	end
	
	if(~isvar('ops') || isempty(ops))
		ops = struct([]);
	end
	
	[y1 P M time hist burn] = part_opt_TRWS_mex(int32(E.G.E-1),double(E.f1),double(E.f2),int32(X0),M,int32(y-1),ops);
else
	fname = varargin{1};
	if(nargin>1)
		M = varargin{2};
	else
		M = [];
	end
	if(nargin>2)
		y = varargin{3};
	else
		y = [];
	end
	if(nargin>3)
		ops = varargin{4};
	else
		ops = struct([]);
	end
	[y1 P M time hist burn] = part_opt_TRWS_mex(fname,M,int32(y-1),ops);
end
x = double(y1)+1;
P = P + 1;
stats.P = double(P);
stats.M = M;
K = size(P,1);
nV = size(P,2);
X = P == repmat([1:K]',1,nV);
%clear part_opt_TRWS_mex;

clear stats;
stats.hist = hist;
stats.burn = burn;
stats.time = time;
stats.LB = hist(1,end);
stats.E = hist(2,end);
stats.tLB = hist(1,:);
stats.tE = hist(2,:);
stats.tt = cumsum(hist(4,:));
stats.tdM = hist(3,:);
end