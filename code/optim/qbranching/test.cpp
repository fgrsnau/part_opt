#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include "optim/qbranching/qbranching.h"


int main(int argc, char *argv[]){
	using namespace exttype;

	if (argc < 1){
		printf("Usage: qbranching edges=edges.cvs edge_pairs=edge_pairs.cvs [option=value]\n");
		printf("No arguments provided");
	} else{
	};
	options ops;
	std::string edges;
	std::string edge_pairs;
	// options are set like this:
	ops["max_CPU"] = 2;
	//ops["max_it"] = 10;
	edges = "testdata/edges.csv";
	edge_pairs = "testdata/edge_pairs.csv";
	//parese input options
	for (int i = 1; i < argc; ++i){
		std::stringstream o(std::string(argv[i]).c_str());
		std::string o_name;
		std::string o_val;
		std::getline(o, o_name, '=');
		std::getline(o, o_val, '=');
		debug::stream << "option " << o_name << " = " << o_val << "\n";
		if (strcmp(o_name.c_str(), "edges")==0){
			edges = o_val;
		} else if (strcmp(o_name.c_str(), "edge_pairs")==0){
			edge_pairs = o_val;
		} else{
			ops[o_name] = atof(o_val.c_str());
		};
	};
	
	dynamic::num_array<int, 2> ee;
	dynamic::num_array<int, 2> hh;
	dynamic::num_array<double, 1> ww;

	for (int loop = 0; loop < 2; ++loop){// first iteration just to count
		int nE = 0;
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
			if (loop == 1){
				ee(0, nE) = u;
				ee(1, nE) = v;
			};
			++nE;
		};
		if (loop == 0){
			ee.resize(mint2(2, nE));
		};
	};

	// Line Graph
	for (int loop = 0; loop < 2; ++loop){
		int nH = 0;
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
				hh(0, nH) = u;
				hh(1, nH) = v;
				ww(nH) = w;
			};
			++nH;
		};
		if (loop == 0){
			hh.resize(mint2(2, nH));
			ww.resize(nH);
		};
	};

	dynamic::num_array<int, 2> B = qbranching(ee, hh, ww, ops);
	// write out solution
	{std::ofstream file(edges + ".sol.csv");
		for (int v = 0; v < B.size()[1]; ++v){
			int u = B(0, v);
			int e = B(1, v);
			file << u << "," << v << "," << e << "\n";
		};
	};
};