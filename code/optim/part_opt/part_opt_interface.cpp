#include "part_opt_interface.h"

#include "dynamic/num_array.h"
#include "energy.h"
#include "part_opt_TRWS.h"

void energy_read_K(int v, int K);

template<typename type>
part_opt_interface<type>::part_opt_interface(){
	energy = new energy_auto<compute_type>();
	alg = 0;
};

template<typename type>
void part_opt_interface<type>::clear_alg(){
	if (alg){
		delete alg;
		alg = 0;
	}
};

template<typename type>
part_opt_interface<type>::~part_opt_interface(){
	clear_alg();
	if (energy){
		delete energy;
	};
};

//loop 0

template<typename type>
void part_opt_interface<type>::energy_read_start(int nV){
	energy->set_nV(nV);
	energy->G._nV = nV;
	energy->K.resize(nV);
	debug::stream << "Converting model\n";
	nE = 0;
};

template<typename type>
void part_opt_interface<type>::energy_read_vertex(int v, int K){
	energy->K[v] = K;
};


template<typename type>
void part_opt_interface<type>::energy_read_edge(int u, int v){
	++nE;
};

template<typename type>
void part_opt_interface<type>::energy_init0(){
	energy->set_nE(nE);
	energy->maxK = energy->K.max().first;
	maxK = energy->maxK;
	nE = 0;
};

//loop 1

template<typename type>
void part_opt_interface<type>::energy_read_f1(int v, int K, double * pval){
	dynamic::num_array<double, 1> f1;
	f1.set_ref(pval, K);
	energy->set_f1(v, f1);
};

template<typename type>
void part_opt_interface<type>::energy_read_f2(int u, int v, int K1, int K2, double * pval){
	dynamic::num_array<double, 2> f2;
	f2.set_ref(pval, exttype::mint2(K1, K2));
	energy->G.E[nE][0] = u;
	energy->G.E[nE][1] = v;
	energy->set_f2(nE, f2);
	++nE;
};

template<typename type>
void part_opt_interface<type>::energy_init1(){
	energy->G.edge_index();
	energy->init();
	energy->report();
	clear_alg();
	alg = new alg_po_trws();
	ops = & alg->ops;
};

template<typename type>
void part_opt_interface<type>::alg_run(){
	alg->set_E(energy);
	debug::PerformanceCounter c1;
	alg->init();
	alg->run_converge();
	iter_init_trws = alg->total_it;
	time_init = c1.time();
	debug::PerformanceCounter c2;
	alg->prove_optimality();
	c2.stop();
	double elim = double(alg->maxelim + alg->nV - alg->nimmovable) / double(alg->maxelim) * 100;
	double time1 = c1.time();
	debug::stream << "INSTANCE " << instance_name << "\n nit: " << alg->total_it << " elim:" << elim << "% / " << time1 << "s" "\n";
	// save statisitcs
	time_total = time1;
	time_po = c2.time();
	iter_po = alg->po_it;
	iter_po_trws = alg->total_it - iter_init_trws;
};

template<typename type>
bool part_opt_interface<type>::is_alive(int v, int k){
	if (!alg){
		throw debug_exception("Must execute Infer first");
	};
	return alg->UU[v][k]; // immovable == alive
};

template class part_opt_interface < float >;
template class part_opt_interface < double >;
template class part_opt_interface < int >;