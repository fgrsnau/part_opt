#include <string>
#include <iostream>
#include <sstream>

#include "dynamic/options.h"
#include "optim/part_opt/part_opt_TRWS.h"
#include "dynamic/dynamic.h"
#include "optim/graph/mgraph.h"
#include "debug/performance.h"

#include "qbranching.h"

dynamic::num_array<int, 2> qbranching(dynamic::num_array<int, 2> & ee, dynamic::num_array<int, 2> & hh, dynamic::num_array<double, 1> & ww, options & ops){
	using exttype::mint2;
	using exttype::mint3;
	using exttype::mint4;

	//Base Graph
	datastruct::mgraph G;
	int nE = ee.size()[1]; // number of edges
	int nV = 0; // number of vertices
	G.E.resize(nE);
	for (int e = 0; e < nE; ++e){
		G.E[e] << ee.subdim<1>(e);
		nV = std::max(nV, G.E[e].max().first + 1);
	};
	G._nV = nV;
	G.edge_index();

	//Line Graph
	datastruct::mgraph H;
	int nH = hh.size()[1]; // number of edges of line graph
	int nL = 0; // number of vertices of line graph = number of edges of G
	H.E.resize(nH);
	for (int h = 0; h < nH; ++h){
		H.E[h] << hh.subdim<1>(h);
		nL = std::max(nL, H.E[h].max().first + 1);
	};
	H._nV = nL;
	H.edge_index();
	
	/*
	//Base Graph
	int nE = 0;
	int nV = 0;	
	for (int loop = 0; loop < 2; ++loop){
		nE = 0;
		nV = 0;
		std::ifstream file(edges);
		if (file.fail()){
			std::cerr << "Error: " << strerror(errno);
			std::cerr << edges.c_str();
		};
		while (file.good()){
			std::string line;
			getline(file, line);
			std::stringstream linest(line);
			int u, v; char c;
			linest >> u >> c >> v;
			if ((linest.rdstate() & std::ifstream::failbit) != 0){
				break;
			};
			//linest >> v;
			if (loop == 1){
				G.E[nE] = mint2(u,v);
			};
			nV = std::max(nV, u + 1);
			nV = std::max(nV, v + 1);
			++nE;
		};
		if (loop == 0){
			G.E.resize(nE);
		} else{
			G._nV = nV;
			G.edge_index();
		};
	};

	// Line Graph
	// model: each vertex selects an incoming edge or the root
	datastruct::mgraph H;
	dynamic::num_array<double, 1> ww; // weight
	int nH = 0;
	int nL = 0; // number of vertices of the line graph = number of edges in G
	for (int loop = 0; loop < 2; ++loop){
		nH = 0;
		nL = 0;
		std::ifstream file(edge_pairs);
		if (file.fail()){
			std::cerr << "Error: " << strerror(errno);
		};
		while (file.good()){
			std::string line;
			getline(file, line);
			std::stringstream linest(line);
			int u, v; char c;
			double w;
			linest >> u >> c >> v >> c >> w;
			if ((linest.rdstate() & std::ifstream::failbit) != 0){
				break;
			};
			if (loop == 1){
				//hh(0, nH) = u;
				//hh(1, nH) = v;
				H.E[nH] = mint2(u, v);
				ww[nH] = w;
			};
			nL = std::max(nL, u + 1);
			nL = std::max(nL, v + 1);
			++nH;
		};
		if (loop == 0){
			//hh.resize(mint2(2, nH));
			H.E.resize(nH);
			ww.resize(nH);
		} else{
			H._nV = nL;
			H.edge_index();
		};
	};
	*/

	// edges of the model are all edges except of the root
	int maxK = 0; // maximum number of labels = max in degree
	int r = -1; // the root
	for (int u = 0; u < nV; ++u){
		maxK = std::max(maxK,G.in[u].size());
		if (G.in[u].size() == 0){
			r = u;
			//happens to be the last one
			assert(r == nV-1);
		};
	};
	// graph without edges from the root
	datastruct::mgraph G1;
	G1._nV = nV - 1;
	G1.E.reserve(nE);
	int g = 0;
	for (int e = 0; e < nE; ++e){
		if (G.E[e][0] != r){//all edges except from the root
			G1.E.push_back(G.E[e]);
		};
	};
	G1.edge_index();
	nV = nV - 1; // set to match G1
	nE = G1.nE(); // set to match G1
	// number of labels = in_degree
	intf K(nV);
	for (int u = 0; u < nV; ++u){
		K[u] = G.in[u].size(); // this includes a label for edge from the "root" = degree in G not in G1
	};
	// pairwise terms
	dynamic::fixed_array1< dynamic::num_array<double, 2> > ff2(nE);
	for (int e = 0; e < nE; ++e){
		int u = G1.E[e][0];
		int v = G1.E[e][1];
		ff2[e].resize(mint2(K[u], K[v]));
		ff2[e] << 0;
	};
	//dynamic::num_array<double, 2> f2(mint2(maxK,maxK));
	for (int h = 0; h < nH; ++h){ // all weighted triplets
		int e1 = H.E[h][0];
		int e2 = H.E[h][1];
		assert(G.E[e1][1] == G.E[e2][0]); // must share a vertex
		// vertices
		int u = G.E[e2][0];
		int v = G.E[e2][1];
		//in_u = G.get_in(u);
		//in_v = G.get_in(v);
		// which number is e1 in the incoming list for u
		int e1i = G.in[u].find(e1);
		assert(e1i >= 0);// found //e1i = find(in_u == e1);
		// which number is e2 in incoming list for v
		//e2i = find(in_v == e2); % cannot be root, in this triplet
		int e2i = G.in[v].find(e2);
		assert(e2i >= 0);// found
		//
		//f2(e1i, e2i, e2) = ww(h);
		ff2[e2](e1i, e2i) = ww(h);
		{	// forbidden assignments : when x(u) = v(and x(v) = u, which is true for this triplet)
			// this forbids the simplest cycle u<->v
			// find v in the incoming list for u
			int e1i = G.find_in(u, v); // e1i = find(G.E(1, in_u) == v);
			ff2[e2](e1i, e2i) = 1000; // forbidden
		};
	};

	// Reformulate without reverse edges
	datastruct::mgraph G2;
	G2._nV = nV;
	G2.E.reserve(nE);
	int nE2 = 0;
	for (int e = 0; e < nE; ++e){
		int u = G1.E[e][0];
		int v = G1.E[e][1];
		if (u < v){ // prefer succesive indexing
			G2.E.push_back(mint2(u, v));
			++nE2;
			int e2 = G1.edge(v, u); // reverse edge
			if (e2 >= 0){ // e2 is the reverese edge u<-v
				for (int i = 0; i < K[u]; ++i){
					for (int j = 0; j < K[v]; ++j){
						ff2[e](i, j) += ff2[e2](j, i);
						//dellist(end + 1) = e2;
					};
				};
			};
		};
	};
	assert(G2.E.size() == nE2);
	for (int e = 0; e < nE2; ++e){
		int u = G2.E[e][0];
		int v = G2.E[e][1];
		assert(u >= 0 && u < nV);
		assert(v >= 0 && v < nV);
	};
	G2.edge_index();
	//
	//Construct energy over graph G2
	energy_auto<float> E;
	E.set_G(G2);
	// number of labels = in_degree
	E.K = K; E.maxK = maxK;
	// no unary terms
	for (int u = 0; u < nV; ++u){
		dynamic::num_array<double, 1> f1(K[u]);
		f1 << 0;
		E.set_f1(u, f1);
		//E.f1[u].resize
		//E.f1[u] << 0;
	};
	for (int e = 0; e < nE2; ++e){
		E.set_f2(e, ff2[e]);
	};
	
	// Energy might want to do some initialization when all data is available
	E.init();
	// Print statistics of the model
	E.report();
	// Create solver
	alg_po_trws<float_v4> alg;
	// Set solver options
	alg.ops << ops;
	// Set energy
	alg.set_E(&E);
	// Init solver
	alg.init();

	debug::PerformanceCounter c1; // start a timer
	// Run TRW-S
	alg.run_converge();
	debug::stream << "TRW-S time:" << c1.time() << "\n";
	// This gives a test labeling y and a starting reparametrization
	// Print labeling
	debug::stream << "Labeling best_x: \n";
	for (int s = 0; s < std::min(nV, 10); ++s){
		debug::stream << alg.best_x[s] << ", ";
	};
	debug::stream << "...\n";
	// Run partial optimality iterations
	alg.prove_optimality();
	double time = c1.time();
	// Percent of eliminated labels
	double elim = alg.elim_labels_percent();
	// Total TRWS iterations
	int nit = alg.total_it;
	debug::stream << " nit: " << nit << " elim:" << elim << "% / " << time << "s" "\n";

	//decode resulting assingment

	dynamic::num_array<int, 2> B(mint2(2, nV));
	for (int u = 0; u < nV; ++u){
		int e1i = alg.best_x(u); // incoming index
		int e1 = G.in[u][e1i];
		B(0, u) = G.E[e1][0];
		B(1, u) = e1;
	};

	return B;

	/*
	% 
		B = zeros(2, nV);
	for u = 1:nV
		e1i = x(u); % incoming index
		in_u = G.get_in(u);
	e1 = in_u(e1i);
	B(1, u) = u;
	B(2, u) = e1;
	end
		%% write solution and log
		dlmwrite([fname '-sol.cvs'], B'-1);
		copyfile('log/output.txt', [fname '-log.txt']);
	% label trees
		G2 = G;
	G2.E = G2.E(:, B(2, :)); % leave only edges from the branching
		G2 = G2.edge_index();
	r = find(G2.in_deg() == 0); % the root, delete it
		G2.V(r) = [];
	G2.E(:, G2.E(1, :) == r) = [];
	G2 = G2.edge_index();
	rr = row(find(G2.in_deg() == 0)); % root nodes
		d_out = G2.out_deg();
	dr = row(d_out(rr)); % print roots degrees
		fid = fopen([fname '-log.txt'], 'a');
	fprintf(fid, '\n%s\n', dd{ di });
	fprintf(fid, 'minimization %3.2f s. (%3.2fs.)\n', toc, stats.time);
	fprintf('TRW-S bound: %d\n', stats.LB);
	fprintf(fid, 'TRW-S lower bound: %d\n', stats.LB);
	fprintf(fid, 'TRW-S solution: %d\n', stats.E);
	fprintf('TRW-S solution: %d\n', stats.E);
	for i = 0:max(dr)
		fprintf('%i trees with root deg %i\n', nnz(dr == i), i);
	fprintf(fid, '%i trees with root deg %i\n', nnz(dr == i), i);
	end
		% label branches
		active = rr;
	label = zeros(1, nV);
	label(rr) = 1:length(rr);
	while ~isempty(active)
		u = active(1);
	active(1) = [];
	for e = G2.get_out(u)
		v = G2.E(2, e); % edge u->v
		label(v) = label(u);
	active(end + 1) = v;
	end
		end
		if (any(label == 0))
			fprintf('cycle detected!\n');
	fprintf(fid, 'cycle detected!\n');
		else
			fprintf('cycle free\n');
	fprintf(fid, 'cycle free\n');
	end
		fclose(fid);
	*/
};
