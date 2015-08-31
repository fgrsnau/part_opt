#include "part_opt_TRWS.h"

#include "exttype/minmax.h"
#include "debug/performance.h"
#include <limits>
#include <algorithm>
//#include <omp.h>

//#include "unroller.h"

using namespace custom_new;

typedef d_type type;
typedef d_vectorizer vectorizer;

//#include "qpbo-v1.3/QPBO.h"
#ifdef WITH_MAXFLOW
#include "graph.h"
using namespace maxflowLib;
#else
#include "maxflow/graph.h"
#endif

using namespace exttype;
using namespace dynamic;


//const double FINF = std::numeric_limits<double>::infinity();
//const int IINF = std::numeric_limits<int>::max();

void err_function(char * s){
	debug::stream << s << "\n";
};

//_________________________alg_po_trws_________________________________

intf alg_po_trws::y0(0);

// construction

alg_po_trws::alg_po_trws(){
	f2_replaced = false;
	parent::ops = &ops;
	burn = 0;
};

alg_po_trws::~alg_po_trws(){ 
	cleanup(); 
};

void alg_po_trws::init1(){
	//total_it = 0;
	po_it = 0;
	nimmovable = 0;
	UU.resize(nV);
	unmark();
	y.resize(nV);
	y << 0;
	parent::set_y(y);
	y1.resize(nV);
};

void alg_po_trws::set_y(intf & y){
	this->y = y;
	parent::set_y(y);
	apply_cut(y, "init");
};

void alg_po_trws::set_P(num_array<int, 2> & P){ // export improving mapping
	this->P.set_ref(P);
};

void alg_po_trws::unmark(){
	for (int v = 0; v < nV; ++v){
		UU[v].reset();
	};
	nimmovable = 0;
};


void alg_po_trws::queue_neib(node_info & v){
	// all neighbors are candidates for node cust
	for (int i = 0; i < v.in.size(); ++i){
		node_info & w = *v.in[i]->tail;
		if (!w.active){
			dee_q.push(&w);
			w.active = true;
		};
	};
	for (int i = 0; i < v.out.size(); ++i){
		node_info & w = *v.out[i]->head;
		if (!w.active){
			dee_q.push(&w);
			w.active = true;
		};
	};
};

int alg_po_trws::mark_immovable(node_info & v, int k, prunner who){
	int s = &v - nodes.begin();
	int nimmovable0 = nimmovable;
	if (ops.all_to_one){
		for (int i = 0; i < v.K; ++i){
			if (!UU[s][i]){
				mark_immovable_label(v, i, who);
			};
		};
	} else{ // subset-to-one
		mark_immovable_label(v, k, who);
	};
	return nimmovable - nimmovable0;
};

void alg_po_trws::rebuild_incremental(node_info & v, int k){
	po_rebuild.resume();
	int s = &v - nodes.begin();
	// rebuild energy unary anew
	po_mask & U_s = UU[s];
	assert(U_s[k] == 1);
	node_info * S = &v;
	if (ops.reduce_immovable && ops.fast_msg){
		if (k == y[s]){
			S->f1[k] = 0;
		} else{
			S->f1[k] = INF(type);
		};
	} else{
		S->f1[k] = 0;
	};
	//rebuild energy pairwise anew
	// outcoming edges
	for (int j = 0; j < v.out.size(); ++j){
		edge_info & ee = *v.out[j];
		int t = ee.head - nodes.begin();
		po_mask & U_t = UU[t];
		if (!U_t[y[t]])continue; // initial rebuild
		ee.f2()->rebuild(U_s, y[s], U_t, y[t]);
		if (ops.reduce_immovable){
			ee.f2()->reduce(U_s, y[s], U_t, y[t]);
		};
	};

	//incoming edges
	for (int j = 0; j < v.in.size(); ++j){
		edge_info & ee = *v.in[j];
		int t = ee.tail - nodes.begin();
		po_mask & U_t = UU[t];
		if (!U_t[y[t]])continue; // initial rebuild
		ee.f2()->rebuild(U_t, y[t], U_s, y[s]);
		if (ops.reduce_immovable){
			ee.f2()->reduce(U_t, y[t], U_s, y[s]);
		};
	};

	best_E = 0;
	po_rebuild.pause();
};

void alg_po_trws::mark_immovable_label(node_info & v, int k, prunner who){
	int s = &v - nodes.begin();
	assert(UU[s][k] == 0);
	if (!UU[s][k]){
		UU[s][k] = 1;
		queue_neib(v);
		++nimmovable;
		// incremental rebuild of the energy
		if (ops.rebuild_incremental){
			rebuild_incremental(v, k);
		};
		// this is message correction for reduction
		if (ops.reduce_immovable){
			int y_s = y[s];
			for (int j = 0; j < v.out.size(); ++j){
				type & m_y = v.out[j]->msg[y_s];
				m_y = std::min(m_y, v.out[j]->msg[k]);
			};
			type & a = v.bw[y_s];
			a = std::min(a, v.bw[k]);
		};
		// save progress for the record
		if (burn){
			(*burn)(k, s, 0) = po_it;
			(*burn)(k, s, 1) = nodes[s].r1[k] - nodes[s].r1[y[s]];
			(*burn)(k, s, 2) = who;
			(*burn)(k, s, 3) = ops.sensetivity;
			(*burn)(k, s, 4) = po_total.time();
			(*burn)(k, s, 5) = parent::total_it;
		};
	};
};

void alg_po_trws::reduce_immovable(){
	for (int e = 0; e < nE; ++e){
		edge_info & ee = edges[e];
		int s = ee.tail - nodes.begin();
		int t = ee.head - nodes.begin();
		po_mask & U_s = UU[s];
		po_mask & U_t = UU[t];
		ee.f2()->reduce(U_s, y[s], U_t, y[t]);
	};//end reduce energy
};

void alg_po_trws::rebuild_energy(){
	po_rebuild.resume();

	//rebuild energy pairwise anew
	for (int e = 0; e < nE; ++e){
		edge_info & ee = edges[e];
		int s = ee.tail - nodes.begin();
		node_info * S = ee.tail;
		int t = ee.head - nodes.begin();
		po_mask & U_s = UU[s];
		po_mask & U_t = UU[t];
		if (!ee.f2())throw debug_exception("pairwise term not initialized");
		ee.f2()->rebuild(U_s, y[s], U_t, y[t]);
	};

	// rebuild energy unary anew
	for (int s = 0; s < nV; ++s){
		po_mask & U_s = UU[s];
		node_info * S = &nodes[s];
		//VALUE_TYPE m1 = FINF;
		for (int k = 0; k < S->K; ++k){
			if (U_s[k]){ // immovable
				if (ops.reduce_immovable && ops.fast_msg){
					if (k == y[s]){
						S->f1[k] = 0;
					} else{
						S->f1[k] = INF(type);
					};
				} else{
					S->f1[k] = 0;
				};
				//m1 = std::min(m1, S->bw[k]);
			} else {
				S->f1[k] = type(E0->f1[s][k] - E0->f1[s][y[s]]) - type(ops.sensetivity);
			};
		};
	};

	if (ops.reduce_immovable){
		reduce_immovable();
	};
	//
	po_rebuild.pause();
	//reset for new energy
	best_E = 0;
};

void alg_po_trws::finish(){
	po_total.stop();
	int elim = 0;
	for (int s = 0; s < nV; ++s){
		for (int k = 0; k < nodes[s].K; ++k){
			if (UU[s][k] == 0){
				++elim;
			};
		};
	};
	if (nimmovable + elim != maxelim + nV){
		throw debug_exception("wrong counts");
	};
	if (P.size()[1] >= nV && P.size()[0] >= maxK){
		debug::stream << "P is right\n";
		P << -1;
		for (int s = 0; s < nV; ++s){
			int y_s = y[s];
			for (int k = 0; k < nodes[s].K; ++k){
				if (UU[s][k]){ //immovable
					P(k, s) = k;
					//X(k, s) = 1; // alive
				} else{
					P(k, s) = y_s;
					//X(k, s) = 0; // eliminated
				};
			};
		};
	};
	debug::stream << "Finished. Total iterations: " << total_it << " Persistency proven: " << double(elim) / maxelim * 100 << "%\n";
	debug::stream << "PO Times:\n Total: " << po_total.time() << "\n Msgs: " << po_msgs.time() << "\n Cuts: " << po_cuts.time() << "\n Node cuts: " << po_node_cuts.time() << "\n Rebuild: " << po_rebuild.time() << "\n";
};

bool alg_po_trws::flat_cuts(){
	intf bbest_x = best_x;
	bool any_cut_succeded = false;
	for (int i = 0; i < maxK; ++i){
		for (int s = 0; s < nV; ++s){
			node_info & v = nodes[s];
			if (i < v.K){
				best_x[s] = i;
			} else{
				best_x[s] = bbest_x[s];
			};
		};
		// run cut_tr
		char s[1024];
		sprintf(s, "flat cut %i", i);
		bool cs = cut_tr(s);
		//if (cs){
		//	debug::stream << "Heuristic cut " << i << " succeeded\n";
		//};
		any_cut_succeded = any_cut_succeded || cs;
	};
	return any_cut_succeded;
};

bool alg_po_trws::cut_tr(char * caller){// make a cut between y and best_x
	// we want to improve best_x by swithching a part of it to y
	// labeling best_x corresponds to the source (label 1)
	// labeling y corresponds to the sink (label 0)
	//
	po_cuts.resume();
	//typedef int tcap;
	//typedef long long tflow;
	typedef type tcap;
	typedef double tflow;
	typedef Graph<tcap, tcap, tflow> tgraph;
	tgraph * g = new tgraph(nV, nE);
	g->add_node(nV);
	//num_array<tgraph::arc*, 1> arcs(nE);
	for (int s = 0; s < nV; ++s){
		tcap g1 = nodes[s].f1[best_x[s]] - nodes[s].f1[y[s]];
		g->add_tweights(s, 2 * g1, 0); // cost g1 if best_x is selected
	};
	int ntrunc = 0;
	for (int e = 0; e < nE; ++e){
		edge_info & ee = edges[e];
		int s = E0->G.E[e][0];
		int t = E0->G.E[e][1];
		tcap E_11 = (*ee.f2())(best_x[s], best_x[t]);
		tcap E_01 = (*ee.f2())(y[s], best_x[t]);
		tcap E_10 = (*ee.f2())(best_x[s], y[t]);
		assert(math::abs((*ee.f2())(y[s], y[t])) < 1e-10);
		tcap E_00 = 0;

		tcap c = E_01 + E_10 - E_11 - E_00;
		if (c < - E0->tolerance){ // not submodular
			throw debug_exception("should be called on the truncated problem already");
		};
		c = std::max(c, tcap(0));
		g->add_edge(s, t, c, c); // cost cap if (x_s,x_t)  is selected
		g->add_tweights(s, 2 * (E_10 - E_00) - c, 0);
		g->add_tweights(t, 2 * (E_01 - E_00) - c, 0);
	};
	tflow flow = g->maxflow(false);
	intf x(nV);
	for (int s = 0; s < nV; ++s){
		bool segment;
		if (ops.weak_po){
			//lower cut, closest to sink -- this way we preserve weak persistencies
			segment = g->what_segment(s, tgraph::SOURCE) == tgraph::SOURCE; // residual connection to the SOURCE tree <=> labeling y is selected
		} else{
			//upper cut, closest to source -- this way we enforce strong persistencies
			segment = g->what_segment(s, tgraph::SINK) == tgraph::SOURCE; // residual connection to the SOURCE tree <=> labeling y is selected
		};

		if (segment){
			x[s] = y[s];
		} else{ // labeling best_x is selected
			x[s] = best_x[s];
		};
	};
	// check what is resulting x compared to best_x -- should be better
	double E_x1 = cost(x);
	debug::stream << "Improved E:" << txt::String::Format("%8.4f", E_x1) << "\n";
#ifdef _DEBUG
	double E_x = cost(best_x);
	if (ops.rebuild_incremental){
		assert(E_x <= 0);
	} else{
		assert(E_x1 <= E_x);
		if (E_x1 <= E_x){
			//debug::stream << "- correct\n";
		} else{
			debug::stream << "- Cut is incorrect\n";
			throw debug_exception("Cut is incorrect");
		};
	};
#endif
	delete g;
	// apply the cut
	po_cuts.pause();
	return apply_cut(x, caller);
};

bool alg_po_trws::apply_cut(intf & x, char * caller){
	int total_mark1 = 0;
	for (int s = 0; s < nV; ++s){
		node_info & v = nodes[s];
		int k = x[s];
		//if (x[s] != y[s] && !UU[s][k]){// mark new cut as immovable
		if (!UU[s][k]){// mark new cut as immovable
			total_mark1 += mark_immovable(v, k, CUT);
		};
	};
	int total_mark2 = mark_dee();
	report_po_iter(total_mark1, total_mark2, caller);
	if (!ops.rebuild_incremental){
		rebuild_energy();
	};
	return total_mark1 > 0;
};

bool alg_po_trws::remove_condition(int s){
	double m = marg(s).first;
	if (ops.weak_po){
		return m < ops.local_min_tol;
	} else {
		return m <= ops.local_min_tol;
	};
};

bool alg_po_trws::stop_condition(){
	for (int s = 0; s < nV; ++s){
		if (remove_condition(s)) return false; // there is something to remove -- do not stop
	};
	return true; //when nothing to remove
};

void alg_po_trws::mark_WTA(){
	// mark all immovable
	int total_mark1 = 0;
	for (int s = 0; s < nV; ++s){
		if (remove_condition(s)){
			int k = marg(s).second; // remove one locally minimal label, that's enough
			total_mark1 += mark_immovable(nodes[s], k, WTA);
		};
	};
	if (total_mark1 == 0 && ops.reduce_immovable)debug_exception("Guaranteed pruning did not happen!");
	int total_mark2 = mark_dee();
	report_po_iter(total_mark1, total_mark2, "WTA");
	if (!ops.rebuild_incremental){
		rebuild_energy();
	};
};

alg_po_trws::t_f2 * alg_po_trws::construct_f2(int e, aallocator * al){
	if (f2_replaced == 0){ // this is not PO phase yet
		return parent::construct_f2(e, al);
	} else{ // this is PO phase, replace pairwise terms with po_reduced
		t_f2 * f2 = 0;
		edge_info & ee = edges[e];
		if (1 && (!ops.fast_msg || (dynamic_cast<term2v_matrix<type, vectorizer>*>(&E0->f2(e))))){ // term is already full matrix or fast_msg disabled
			// use full matrix for po
			return al->allocate<term2v_matrix_po<type, vectorizer> >(E0->f2(e));
		};
		// fast_msg enabled and term is not a full matrix
		if (dynamic_cast<term2v_potts<type, vectorizer>*>(&E0->f2(e))){
			return al->allocate<term2v_po_reduced<term2v_potts<type, vectorizer> > >(dynamic_cast<term2v_potts<type, vectorizer> & >(E0->f2(e)));
		};
		if (dynamic_cast<term2v_tlinear<type, vectorizer>*>(&E0->f2(e))){
			return al->allocate<term2v_po_reduced<term2v_tlinear<type, vectorizer> > >(*dynamic_cast<term2v_tlinear<type, vectorizer>*>(&E0->f2(e)));
		};
		if (dynamic_cast<term2v_diff<type, vectorizer>*>(&E0->f2(e))){
			return al->allocate<term2v_po_reduced<term2v_diff<type, vectorizer> > >(*dynamic_cast<term2v_diff<type, vectorizer>*>(&E0->f2(e)));
		};
		if (dynamic_cast<term2v_tquadratic<type, vectorizer>*>(&E0->f2(e))){
			return al->allocate< term2v_po_reduced<term2v_tquadratic<type, vectorizer> > >(*dynamic_cast<term2v_tquadratic<type, vectorizer>*>(&E0->f2(e)));
		};
		if (dynamic_cast<term2v_matrix<type, vectorizer>*>(&E0->f2(e))){
			return al->allocate<term2v_po_reduced<term2v_matrix<type, vectorizer> > >(*dynamic_cast<term2v_matrix<type, vectorizer>*>(&E0->f2(e)));
		};
		throw debug_exception("Problem: cannot find the true type of the pairwise term");
	};
};

void alg_po_trws::init_aux(){
	// need to recreate the core with substituted pairwise terms
	// first save messages
	dynamic::fixed_array1<tvect> M(nE);//backward messages
	for (int e = 0; e < nE; ++e){
		int s = E0->G.E[e][0]; //tail node
		int K = nodes[s].K;
		M[e].resize(K);
		M[e] << edges[e].msg;
	};
	// destroy core
	destroy_core();
	// do something to substitute new terms
	f2_replaced = true;
	// recreate core with new pairwise terms allocated in the core
	init_core();
	parent::set_y(y);
	// set messsages
	for (int e = 0; e < nE; ++e){
		edges[e].msg << M[e];
	};
	init_iteration(); // restore bw sums
};
void alg_po_trws::cleanup(){
};

void alg_po_trws::set_E(tenergy * E){
	this->E0 = E;
	nV = E->G.nV();
	nE = E->G.nE();
	maxK = E->maxK;
};

void alg_po_trws::init(){
	parent::init(E0);
	maxelim = E0->K.sum() - nV;
	init1();
	ops.print();
};

bool alg_po_trws::check_zero_gap(){
	if (exit_reason == zero_gap){ // TRW-S found optimal solution, nothing to do
		debug::stream << "TRW-S found optimal solution\n";
	};
	if (stop_condition()){
		if (ops.weak_po){
			debug::stream << "Map is weakly improving\n";
		} else{
			debug::stream << "Map is strictly improving\n";
		};
		finish();
		return true;
	} else{
		if (ops.weak_po){
			debug::stream << "Map is not yet weakly improving\n";
		} else{
			debug::stream << "Map is not yet strictly improving\n";
		};
	};
	return false;
};

std::pair<double,int> alg_po_trws::marg(int s){
	std::pair<double, int> r;
	r.second = - 1;
	double m = INF(double);
	po_mask & U_s = UU[s];
	for (int i = 0; i < nodes[s].K; ++i){ 
		if (U_s[i]) continue; // this one is excluded (also with y)
		double v = nodes[s].r1[i];
		if (v < m){
			m = v;
			r.second = i;
		};
	};
	r.first =  m - nodes[s].r1[y[s]]; // minimum of all alive - value at y
	return r;
};

void alg_po_trws::prove_optimality(intf & y0){
	// bla
	//int maxit = ops.po_max_it;
	//debug::stream << "PO Options:\n";
	//debug::stream << "conv_tol: " << ops.conv_tol << "\t    max_it: " << maxit << "\n";
	//debug::stream << "  reduce: " << ops.reduce_immovable << "\t pixel_cut: " << ops.use_pixel_cut << "\n";
	//debug::stream << "    cuts: " << ops.use_cut << "\t  fast_msg: " << ops.fast_msg << "\n";
	//debug_LB = true;
	po_total.start();
	po_cuts.start(); po_cuts.pause();
	po_msgs.start(); po_msgs.pause();
	po_rebuild.start(); po_rebuild.pause();
	po_node_cuts.start(); po_node_cuts.pause();
	if (!y0.empty()){
		y = y0;
		parent::set_y(y);
	} else{
		//assume we are running on a converged state
		y = best_x;
		parent::set_y(y);
		if (check_zero_gap()) return;
		// once y is fixed, there might be single-pixel cuts, if allowed improve them
		icm();
	};

	init_aux();

	// initial aux energy
	double max_sens = -INF(double);
	for (int s = 0; s < nV; ++s){
		node_info & v = nodes[s];
		int k = y[s];
		//mark_immovable_label(v, k, INIT);
		for (int i = 0; i < v.K; ++i){
			double d = (v.r1[k] - v.r1[i]) / v.core->coeff; // this should guarantee 100 acceptance
			max_sens = std::max(max_sens, d);
		};
	};
	max_sens *= 2;
	/*
	for (int e = 0; e < nE; ++e){
		edge_info & ee = edges[e];
		node_info * s = ee.tail;
		node_info * t = ee.head;
		int k = y[s - nodes.begin()];
		for (int i = 0; i < t->K; ++i){
			max_sens = std::max(max_sens, double(ee.msg[k] - ee.msg[i]));
		};
	};
	*/
	// set initial sensetivity to allow super weak persistency
	if (ops.n_sensetivity > 1){
		ops.sensetivity = -max_sens;
	};

	// init unaries 
	for (int s = 0; s < nV; ++s){
		node_info & v = nodes[s];
		for (int i = 0; i < v.K; ++i){
			v.f1[i] = type(E0->f1[s][i] - E0->f1[s][y[s]]) - type(ops.sensetivity);
		};
	};

	apply_cut(y, "init");

	po_it = 1;
	double target_E = -INF(double);
	if (ops.use_cut){
		target_E = 0;
	};
	// switch to PO options
	ops.max_it = ops.po_max_it;
	ops.it_batch = ops.po_it_batch;
	// debug check
	check_get_col();

	for (int si = 0; si<ops.n_sensetivity; ++si, ops.sensetivity = ops.sensetivity + max_sens / (ops.n_sensetivity - 1)){
		debug::stream << "PO_sensetivity:" << ops.sensetivity << "\n";
		for (;; ++po_it){
			debug::stream << "PO_it:" << po_it << "\n";
			// run more TRW-S iterations
			po_msgs.resume();
			run_converge(target_E);
			po_msgs.pause();
			//
			// do some pruning
			if (exit_reason == zero_gap){ // found pruning cut
				// still check to keep more labels
				if (ops.use_cut){
					if (cut_tr("cut")) continue; // has eliminated somethis
				};
			};
			if (exit_reason == target_energy){ // found pruning cut
				if (ops.use_cut){
					cut_tr("cut"); // guaranteed to succeed, cut, dee and rebuild
				};
				continue;
			};
			if (check_zero_gap()) break;
			// if got to here, then exit reason is convergence, iteration limit or there are ties in the optimal solution
			mark_WTA();
		};//end outer iteration
		if (!ops.rebuild_incremental){
			rebuild_energy();
		};
	}; // end sensetivity loop
};

void alg_po_trws::report_po_iter(int total_mark1, int total_mark2, std::string caller){
	debug::stream << "Marking immovable (" << caller << ") " << total_mark1 << "+" << total_mark2;
	if (dee_mark_strong){
		debug::stream << " (" << dee_mark_strong << " non-unique)";
	};
	debug::stream << " labels (" << double(maxelim + nV - nimmovable) / double(maxelim) * 100.0 << "% remains)\n";
};

double alg_po_trws::dee_delta(node_info & v, int k1){
	int s = &v - nodes.begin();
	int y_s = y[s];
	double Delta = E0->f1[s][k1] - E0->f1[s][y_s] - ops.sensetivity;
	for (int j = 0; j < v.out.size(); ++j){
		edge_info * ee = v.out[j];
		int e = ee - edges.begin();
		node_info * w = ee->head;
		int t = w - nodes.begin();
		type a = INF(type); // std::numeric_limits<VALUE_TYPE>::max();
		for (int k2 = 0; k2 < w->K; ++k2){
			if (!UU[t][k2])continue; // minimize only over immovable
			type V = E0->f2(e)(k1, k2) - E0->f2(e)(y_s, k2);
			//assert(V == floor(V));
			//int v = int(V);
			if (V < a){
				a = V;
				//y1[t] = k2;
			};
		};
		Delta += a;
	};
	for (int j = 0; j < v.in.size(); ++j){
		edge_info * ee = v.in[j];
		int e = ee - edges.begin();
		node_info * w = ee->tail;
		int t = w - nodes.begin();
		type a = INF(type);
		for (int k2 = 0; k2 < w->K; ++k2){
			if (!UU[t][k2])continue; // minimize only over immovable
			type V = E0->f2(e)(k2, k1) - E0->f2(e)(k2, y_s);
			//assert(V == floor(V));
			//int v = int(V);
			if (V < a){
				a = V;
				//y1[t] = k2;
			};
		};
		Delta += a;
	};
	return Delta;
};

double alg_po_trws::icm_delta(node_info & v, int k1){
	int s = &v - nodes.begin();
	int y_s = y[s];
	double Delta = E0->f1[s][k1] - E0->f1[s][y_s];
	for (int j = 0; j < v.out.size(); ++j){
		edge_info * ee = v.out[j];
		int e = ee - edges.begin();
		node_info * w = ee->head;
		int t = w - nodes.begin();
		type V = E0->f2(e)(k1, y[t]) - E0->f2(e)(y_s, y[t]);
		Delta += V;
	};
	for (int j = 0; j < v.in.size(); ++j){
		edge_info * ee = v.in[j];
		int e = ee - edges.begin();
		node_info * w = ee->tail;
		int t = w - nodes.begin();
		type V = E0->f2(e)(y[t], k1) - E0->f2(e)(y[t], y_s);
		Delta += V;
	};
	return Delta;
};

void alg_po_trws::icm(){
	for (int s = 0; s < nV; ++s){
		node_info & w = nodes[s];
		if (!w.active){
			dee_q.push(&w);
			w.active = true;
		};
	};
	while (!dee_q.empty()){
		node_info & v = *dee_q.front();
		v.active = false;
		dee_q.pop();
		int s = &v - nodes.begin();
		for (int k1 = 0; k1 < v.K; ++k1){
			double Delta = icm_delta(v, k1);
			//break;
			if (Delta < 0){// exist immovable labeling of neighbors such that mapping label k is not strictly improving
				y[s] = k1;
				queue_neib(v);
				//queue self back as well
				if (!v.active){
					dee_q.push(&v);
					v.active = true;
				};
				//debug::stream << "ICM improved E: " << txt::String::Format("%f", cost(y)) << "\n";
			};
		};
	};
	best_E = cost(y);
	best_x = y;
	parent::set_y(y);
	debug::stream << "ICM improved E: " << txt::String::Format("%f", best_E) << "\n";
};

bool alg_po_trws::is_alive(int s, int k){
	return UU[s][k];
};

int alg_po_trws::test_label(int s){
	return y[s];
};

double alg_po_trws::elim_labels_percent(){
	return double(maxelim + nV - nimmovable) / double(maxelim) * 100;
};

int alg_po_trws::mark_dee(){
	// mark more immovable DEE style
	dee_mark_strong = 0; // count marks
	if (!ops.use_pixel_cut){
		return 0;
	};
	po_node_cuts.resume();
	int total_mark2 = 0;
	while (!dee_q.empty()){
		node_info & v = *dee_q.front();
		v.active = false;
		dee_q.pop();
		//y1 << y;
		int s = &v - nodes.begin();
		for (int k1 = 0; k1 < v.K; ++k1){
			if (UU[s][k1])continue;
			//y1[s] = k1;
			double Delta = dee_delta(v, k1);
			//break;
			/*
			if (Delta <= 0){
			double Ey1 = E0->cost(y1);
			y1[s] = y[s]; // switch
			double Ey2 = E0->cost(y1);
			if (Ey1 - Ey2 > 0){
			debug_exception("something is wrong with pixel prunning\n");
			};
			};
			*/
			//
			if (Delta <= 0 && Delta >= ops.local_min_tol){
				//continue;
				if (ops.weak_po)continue;
				++dee_mark_strong;
			};
			if (Delta <= 0){// exist immovable labeling of neighbors such that mapping label k is not strictly improving
				//v.U[k1] = 1;
				total_mark2 += mark_immovable(v, k1, PXCUT);
			};
		};
	};
	po_node_cuts.pause();
	return total_mark2;

	/*
do{
nreduce = 0;
//break;
for (int s = 0; s < nV; ++s){
node_info & v = nodes[s];
double Delta = 0;
int y_s = y[s];
for (int k1 = 0; k1 < v.K; ++k1){
if (UU[s][k1])continue;
double Delta = E0->f1[s][k1] - E0->f1[s][y_s];
for (int j = 0; j < v.out.size(); ++j){
edge_info * ee = v.out[j];
int e = ee - edges.begin();
node_info * w = ee->head;
int t = w - nodes.begin();
int a = IINF;
for (int k2 = 0; k2 < w->K; ++k2){
if (!UU[t][k2])continue; // minimize only over immovable
VALUE_TYPE V = E0->f2(e)(k1, k2) - E0->f2(e)(y_s, k2);
assert(V == floor(V));
int v = int(V);
a = std::min(a, v);
};
Delta += a;
};
for (int j = 0; j < v.in.size(); ++j){
edge_info * ee = v.in[j];
int e = ee - edges.begin();
node_info * w = ee->tail;
int t = w - nodes.begin();
int a = IINF;
for (int k2 = 0; k2 < w->K; ++k2){
if (!UU[t][k2])continue; // minimize only over immovable
VALUE_TYPE V = E0->f2(e)(k2, k1) - E0->f2(e)(k2, y_s);
assert(V == floor(V));
int v = int(V);
a = std::min(a, v);
};
Delta += a;
};
if (Delta <= 0){// exist immovable labeling of neighbors such that mapping label k is not strictly improving
//v.U[k1] = 1;
mark_immovable(v, k1);
++nreduce;
};
};
};
total_mark2 += nreduce;
} while (nreduce > 0); // end reduce

};
po_node_cuts.pause();
return total_mark2;
*/
};