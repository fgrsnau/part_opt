#include <string>
#include <iostream>
#include <sstream>

#include "dynamic/options.h"
#include "part_opt_TRWS.h"
#include "dynamic/dynamic.h"
#include "optim/graph/mgraph.h"
#include "debug/performance.h"


void test_rand(int rand_inst, int ptype, options & ops, int & nit, double & elim, double & time){
	using exttype::mint2;
	using exttype::mint3;
	using exttype::mint4;

	datastruct::mgraph G;
	int M = 320;
	int N = 240;
	int K = 12;
#ifdef _DEBUG
	M = 20;
	N = 20;
	K = 9;
#endif
	// Construct grid graph
	G.create_grid<2, 0>(mint2(M, N));
	G.edge_index();
	int nV = G.nV();
	int nE = G.nE();

	//Construct energy over graph G
	energy_auto<d_type> f;
	f.set_G(G);
	//Allow K labels in every vertex
	f.K << K;
	//the maximum number of labels over all vertices
	f.maxK = K;

	// Initialize random seed, for a reproducible test
	srand(rand_inst);

	//generate random unaries
	for (int s = 0; s < nV; ++s){//loop over vertices
		dynamic::num_array<double, 1> f1(K);
		for (int k = 0; k < K; ++k){
			f1(k) = double(rand() % 10000) / 100;
		};
		//set to energy
		f.set_f1(s, f1);
	};

	for (int e = 0; e < nE; ++e){// loop over edges
		dynamic::num_array<double, 2> f2(mint2(K, K));
		// there are several types of models that can be tried
		if (ptype == 3){// generate random truncated quadratic model
			double gamma = (rand() % 6000) / 100 + 1;
			double th = ((rand() % 4000 + 3000)*gamma) / 100;
			for (int k1 = 0; k1 < K; ++k1){
				for (int k2 = 0; k2 < K; ++k2){
					double v = std::min(math::sqr(k1 - k2)*gamma, th);
					f2(k1, k2) = v;
				};
			};
		};
		if (ptype == 2){// generate random truncated linear model
			int gamma = rand() % 20 + 1;
			int th = (rand() % 2 + 2)*gamma;
			for (int k1 = 0; k1 < K; ++k1){
				for (int k2 = 0; k2 < K; ++k2){
					f2(k1, k2) = std::min(std::abs(k1 - k2)*gamma, th);
				};
			};
		};
		if (ptype == 1){// generate random potts model
			double  gamma = double(rand() % 5500) / 100;
			f2 << gamma;
			for (int k = 0; k < K; ++k){
				f2(k, k) = 0;
			};
		};
		if (ptype == 0){// generate random full model
			f2 << 0;
			for (int k1 = 0; k1 < K; ++k1){
				for (int k2 = 0; k2 < K; ++k2){
					f2(k1, k2) = double(rand() % 5000) / 100;
				};
			};
		};
		// set to energy, if matrix f2 is of a particular type, it will be recognized
		f.set_f2(e, f2);
	};
	// Energy might want to do some initialization when all data is available
	f.init();
	// Print statistics of the model
	f.report();
	// Create solver
	alg_po_trws alg;
	// Set solver options
	alg.ops << ops;
	// Set energy
	alg.set_E(&f);
	// Init solver
	alg.init();
	
	debug::PerformanceCounter c1; // start a timer
	// Run TRW-S
	alg.run_converge();
	debug::stream << "TRW-S time:" << c1.time() << "\n";
	// This gives a test labeling y and a starting reparametrization
	// Print labeling
	debug::stream << "Labeling best_x: \n";
	for (int s = 0; s < std::min(nV, 10);++s){
		debug::stream << alg.best_x[s] << ", ";
	};
	debug::stream << "...\n";
	// Run partial optimality iterations
	alg.prove_optimality();
	time = c1.time();
	// Find out which labels were eliminated
	debug::stream << "Alive Labels mask: \n";
	for (int k = 0; k < K; ++k){
		debug::stream << " ";
		for (int s = 0; s < std::min(nV, 10); ++s){
			debug::stream << alg.is_alive(s,k) << ", ";// UU is the set of immovable labes
		};
		debug::stream << "...\n";
	};
	//labeling again
	debug::stream << "Labeling y: \n"; //may differ by the ICM improvement
	for (int s = 0; s < std::min(nV, 10); ++s){
		debug::stream <<"(" << alg.test_label(s) << ")";
	};
	debug::stream << "...\n";
	// Percent of eliminated labels
	elim = alg.elim_labels_percent();
	// Total TRWS iterations
	nit = alg.total_it;
	debug::stream << "TYPE:" << ptype << " TEST INSTANCE " << rand_inst << " nit: " << nit << " elim:" << elim << "% / " << time << "s" "\n";
};

int main(int argc, char *argv[]){
	if (argc < 1){
		printf("Usage: test_random [option=value]\n");
		printf("No arguments provided");
	} else{
	};
	options ops;
	// options are set like this:
	ops["max_CPU"] = 2;
	//ops["max_it"] = 10;
	//parese input options
	for (int i = 1; i < argc; ++i){
		std::stringstream o(std::string(argv[i]).c_str());
		std::string o_name;
		std::string o_val;
		std::getline(o, o_name, '=');
		std::getline(o, o_val, '=');
		debug::stream << "option " << o_name << " = " << o_val << "\n";
		ops[o_name] = atof(o_val.c_str());
	};
	//random testsL
	for (int inst = 1; inst <= 1; ++inst){
		for (int ptype = 0; ptype < 4; ++ptype){
			int nit1; double elim1; double time1;
			test_rand(inst, ptype, ops, nit1, elim1, time1);
			debug::stream << "press a key\n";
			std::cin.get();
		};
	};
	//std::cin.get();
};