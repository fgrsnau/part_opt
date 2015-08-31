#include "trws_machine.h"
#include <limits>
#include <stdlib.h>


#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#else
#define _alloca alloca
#endif

//#if defined(USE_OPENMP) || defined(_OPENMP)
//	#include <omp.h>
//#endif


#define PREFETCH(p) _mm_prefetch((char*)(p), _MM_HINT_T1);
#define PREFETCH0(p) _mm_prefetch((char*)(p), _MM_HINT_T0);
#define CL1 64

template<class type, class vectorizer>
template<bool dir_fw, bool measure> void trws_machine<type, vectorizer>::line_update(const node_info_core * v, const node_info_core * v_end, line_compute_struct<type> & temp){
	const int V = sizeof(vectorizer) / sizeof(type);
	int maxKa = v_align<type, vectorizer>(maxK);
	vectorizer * const t1 = (vectorizer*)_alloca(sizeof(type)*maxKa);
	vectorizer * const t2 = measure ? (vectorizer*)_alloca(sizeof(type)*maxKa) : 0;
	vectorizer * const t3 = measure ? (vectorizer*)_alloca(sizeof(type)*maxKa) : 0;
	vectorizer * const t4 = measure ? (vectorizer*)_alloca(sizeof(type)*maxKa) : 0;
	vectorizer * const temp_r1 = measure ? 0 : (vectorizer*)_alloca(sizeof(type)*maxKa);
	tvect t3r;
	//
	int v_dist = 0;
	if (v->next){
		v_dist = (const char *)v->next - (const char *)v;
	};
	//const int c_dist = v_dist;
	const int chunk_dist = (v_dist+CL1-1) / CL1;
	// go few nodes ahead for prefetch
	const node_info_core * w = v;
	for (int i = 0; i < 2; ++i){ // 2 nodes ahead
		if (w->next){
			w = w->next;
		} else{
			break;
		};
	};
	//
	for (; v != v_end; v = v->next){
		assert(v != 0);
		const int vK = v->vK;
		//vectorizer * const pr = v->r1();

		// prefetch
		if (w != 0){ // could go past the end as well, to the next line of nodes
			const char *p = (const char *)(w->next);
			for (int i = 0; i < chunk_dist; ++i){
				PREFETCH(p + i * CL1);
			};
			if (!dir_fw){// backward pass: outcoming mesasges are cold
				for (int j = 0; j < w->template dir_n_out<dir_fw>(); ++j){
					const edge_info_core & e = w->template dir_out<dir_fw>(j);
					PREFETCH(e.msg);
					PREFETCH((char *)e.msg + CL1);
					PREFETCH(e.f2);
					PREFETCH((char *)e.f2 + CL1);
				};
			} else{
			};
			w = w->next;
		};

		// min-marginals will be temporary or stored for a measuring iteration
		vectorizer * const r1 = measure ? v->r1() : temp_r1;
		// vectorizer * const r1 = v->r1();

		if (measure){ 
			//save current min-marginals in t1
			//save bw+f1 in t3, f1 in t4
			for (int i = 0; i < vK; ++i){
				//t1[i] = pr[i];
				t1[i] = v->r1()[i];
				t4[i] = v->f1()[i];
				if (dir_fw){
					t3[i] = v->bw()[i] + v->f1()[i]; // in forward pass must add f1
				} else{
					t3[i] = v->bw()[i]; // in backward pass has f1 already included
				};
			};
#ifdef _DEBUG
			if (!dir_fw){
				for (int i = v->node->K; i < vK*V; ++i){
					assert(((type*)t3)[i] == INF(type));
				};
			};
#endif
		};

		// compute the sum of all incoming messages (forward and backward) and the unary term
		vectorizer norm(INF(type));
		{
			const vectorizer coeff(v->coeff);
			const edge_info_core * const in_beg = v->template dir_in_begin<dir_fw>();
			const edge_info_core * const in_end = v->template dir_in_begin<dir_fw>() + v->template dir_n_in<dir_fw>();
			const vectorizer * pf1 = v->f1();
			vectorizer * const pbw = v->bw();
			for (int i = 0; i < vK; ++i){
				vectorizer a;
				if (dir_fw){
					a = pf1[i];
				} else{
					a = 0;
				};
				for (const edge_info_core * e = in_beg; e != in_end; ++e){
						a += e->msg[i];
						assert(v_first(a) > -INF(type));
						PREFETCH0((const char *)&(e->msg[i]) + sizeof(vectorizer));
				};

				vectorizer b = a + pbw[i];
				b *= coeff;
				v_min(norm, b);
				pbw[i] = a; //after forward pass holds: f1+incoming, after bw pass: only incoming
				r1[i] = b;
			};
			v_min(norm); // min over entries of norm
#ifdef _DEBUG
			if (dir_fw){
				for (int i = v->node->K; i < vK*V; ++i){
					assert(((type*)pbw)[i] == INF(type));
				};
			};
#endif
		};

		if (measure){
			double dLB = double(v_first(norm)) / v->coeff;
			assert(dLB == dLB); //check nan
			temp.LB += dLB;
		};
	
		if (measure){ // compute LB, renormalize, etc., only once in a while; not optimized
			// compute min
			// estimate a better labeling serial rounding
			// t3 holds sum of outcoming + f1 
			// t4 is f1 only
			for (int j = 0; j < v->template dir_n_in<dir_fw>(); ++j){
				const edge_info_core & e = v->template dir_in<dir_fw>(j);
				{
					// normalize incoming messages, helps for numerical stability
					vectorizer m(INF(type));
					for (int i = 0; i < vK; ++i){
						v_min(m, e.msg[i]);
					};
					v_min(m);
					for (int i = 0; i < vK; ++i){
						e.msg[i] -= m;
					};
				};
				// get tail node, get vector
				const node_info_core * tc = e.template dir_tail<dir_fw>();
				//const int vKt = tc->vK;
				int t_x = tc->x;
				assert(t_x >= 0 && t_x < tc->node->K);
				assert(tc->node->f1[t_x] <  INF(type));
				if (dir_fw){
					tvect t2r;
					t2r.set_ref(t2, vK*V);
					e.f2->get_row(t_x, t2r);
				} else{
					e.f2->get_col_e(vK, t_x, t2);
				};
				// aggregate in t4 and t3
				for (int i = 0; i < vK; ++i){
					t4[i] += t2[i];
					t3[i] += t2[i];
				};
			};

#ifdef _DEBUG
			if (!dir_fw){
				for (int i = v->node->K; i < vK*V; ++i){
					assert(((type*)t3)[i] == INF(type));
				};
			};
#endif
			 // t3 is now sum of outcoming + f1 + incoming(x), minimize to get x
			//node_info & node = *v->node;
			t3r.set_ref(t3, vK * V);
			//int y_s = y[int(v->node - nodes.begin())];
			const int y_s = v->y;
			{// find best v->x , prefer y
				std::pair<type, int> a = t3r.min(); //temp.t3.min();
				(const_cast<node_info_core*>(v))->x = a.second;
				// selected new x
				// add to energy
				temp.E += ((const type*)t4)[a.second]; // f1(x) + incoming edges (x)
				assert(v->x >= 0 && v->x < v->node->K);
				assert(v->node->f1[v->x] < INF(type));
			};
			// margin
			type a_first = v_first(norm);
			//type mrg = a_first - ((type*)pr)[y_s];// / v->coeff;
			type mrg = a_first - ((type*)r1)[y_s];// / v->coeff;
			temp.margin = std::min(temp.margin, mrg);
			// dPhi (change of messages)
			vectorizer norm_t1(INF(type));
			for (int i = 0; i < vK; ++i){
				v_min(norm_t1, t1[i]);
			};
			v_min(norm_t1);
			norm_t1 -= norm;
			vectorizer vm(INF(type));
			for (int i = 0; i < vK; ++i){
				vectorizer a = r1[i];
				//vectorizer a = v->r1()[i];
				//v_min(a, vectorizer(1e30));// de-INF a
				vectorizer b = t1[i];
				//v_min(b, vectorizer(1e30));// de-INF b
				a -= b;
				a += norm_t1;
				v_min(a,-a); // -abs value
				v_min(a,vm); // min(-abs)
				vm = a;
			};
			v_min(vm);
			assert(!(v_first(vm) != v_first(vm))); // check not NAN
			temp.dPhi = std::max(temp.dPhi, -v_first(vm) / v->coeff);
		};
		assert(v_first(norm) > -INF(type));
		// Push outcoming messages
		{
			const edge_info_core * const out_beg = v->template dir_out_begin<dir_fw>();
			const edge_info_core * const out_end = v->template dir_out_begin<dir_fw>() + v->template dir_n_out<dir_fw>();
			for (const edge_info_core * e = out_beg; e != out_end; ++e){
				//vectorizer * const pr = v->r1();
					{
						const vectorizer * pm = e->msg;
						for (int i = 0; i < vK; ++i){
							t1[i] = r1[i] - pm[i] - norm;
							//t1[i] = pr[i] - pm[i] - norm;
							assert(v_first(t1[i]) > -INF(type));
							//if (!dir_fw){
							//PREFETCH0((const char *)&pm[i] + c_dist);
							//PREFETCH0((const char *)&pm[i] + c_dist);
							//};
						};
					};
				vectorizer * pm = e->msg;
				if (dir_fw){
					e->f2->min_sum_e(vK, t1, pm);
				} else{
					e->f2->min_sum_et(vK, t1, pm);
				};
			};
		};
	};
};

template<class type, class vectorizer>
template<bool dir_fw, bool measure> void trws_machine<type, vectorizer>::run_dir(){
	if (measure){
		LB = 0;
		dPhi = 0;
		margin = INF(type);
	};

	int NT = get_num_chunks();

	dynamic::fixed_array1<line_compute_struct<type> > temps(NT);
	dynamic::fixed_array1<line_compute_struct<double> > temps_block(NT);
	for (int t = 0; t < NT; ++t){
		temps[t].init(maxK);
		temps_block[t].init(0);
	};
#pragma omp parallel  num_threads(NT) //if(SG.nn.size() > 4)
	{
		int l;
		if (dir_fw){
			// forward loop: l=0; l< SG.nn.size(); ++ l
			l = 0;
		} else{
			// backward loop: l = SG.nn.size()-1; l>=0; --l
			l = SG.nn.size() - 1;
		};
		int th_id = 0;
		// layer loop (forward or backward)
		for (;;){
#pragma omp for schedule(static,1) nowait// one thread per loop variable = NT threads
			for (int th_id = 0; th_id < NT; ++th_id){
				chunk_iterator ci(this);
				ci.th_id = th_id;
				slice_iterator si(ci);
				si.l = l;

				node_info * S = si.first();
				node_info * T = si.last();
				if (S){
					assert(T);
					line_update<dir_fw, measure>(S->core, T->core->next, temps[th_id]);
					if (measure){
						temps_block[th_id].aggregate(temps[th_id]);
					};
				};
			};
#pragma omp barrier // advance to the next layer, only after the previous layer is fully complete
			{
				//std::cout << "\n";
				if (dir_fw){
					if (++l == SG.nn.size())break;
				} else{
					if (--l < 0)break;
				};
			};
		};
	};
	// gather LB
	if (measure){
		for (int t = 1; t < NT; ++t){
			temps_block[0].aggregate(temps_block[t]);
		};
		current_E = temps_block[0].E;
		LB = temps_block[0].LB;
		dPhi = temps_block[0].dPhi;
		margin = temps_block[0].margin;
	};
};

template<class type, class vectorizer>
template<bool measure> void trws_machine<type, vectorizer>::run_forward_backward(){
	run_dir<true, false>();
	run_dir<false, measure>(); // only backward pass measures
	int z = 0;
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::set_y(const intf & y){
	for (int v = 0; v < nV; ++v){
		nodes[v].core->y = y[v];
	};
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::init_iteration(){//! do a backward pass to restore partial sums in bw, assume given messages from a backward pass
	for (int l = SG.nn.size() - 1; l >= 0; --l){//dummy backward pass
		int nl = SG.nn[l];
		for (int i = 0; i < nl; ++i){
			int s = SG.levels[SG.level_start[l]+i];
			node_info & v = nodes[s];
			v.bw << 0;
			for (int j = 0; j < v.out.size(); ++j){
				edge_info * e = v.out[j];
				v.bw += e->msg;
			};
		};
	};
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::renormalize_to_y(){//! normalization such that message is zero at y
	for (int l = SG.nn.size() - 1; l >= 0; --l){//backward pass
		int nl = SG.nn[l];
		for (int i = 0; i < nl; ++i){
			int s = SG.levels[SG.level_start[l] + i];
			int y_s = nodes[s].core->y;
			node_info & v = nodes[s];
			for (int j = 0; j < v.out.size(); ++j){
				edge_info * e = v.out[j];
				e->msg -= e->msg[y_s];
			};
		};
	};
};


template<class type, class vectorizer>
void trws_machine<type, vectorizer>::run(int niter){
	//data_total = 0;
	// run niter-1 without measure
	for (int it = 0; it < niter - 1; ++it){
		run_forward_backward<false>();
	};
	// run once with measure
	run_forward_backward<true>();
	//read off labeling
	for (int v = 0; v < nodes.size(); ++v){
		current_x[v] = nodes[v].core->x;
	};
	total_it += niter;
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::run_converge(double target_E){
	int maxit = ops->max_it;
	exit_reason = error;
	best_E = INF(type);
	int Max_it = int(ceil(double(maxit) / ops->it_batch));
	hist.resize(mint2(4, Max_it));
	double total_t1 = 0;
	int it = 0;
	//init_iteration();
	double LB0 = -INF(double);
	for (; it < Max_it; ++it){
		debug::PerformanceCounter c1;
		run(ops->it_batch);
#ifdef _DEBUG
		double current_E_ref = cost(current_x); // optimize
		assert(std::abs(current_E - current_E_ref) < 1e-4);
#endif
		if (current_E < best_E){
			best_E = current_E;
			best_x = current_x;
			//int bi = sols.size();
			//sols.resize(sols.size() + 1);
			//sols[bi].resize(nV);
			//sols[bi] << best_x;
		};
		double t1 = c1.time();
		total_t1 += t1;

		hist[mint2(0, it)] = LB;
		hist[mint2(1, it)] = current_E;
		hist[mint2(2, it)] = dPhi;
		hist[mint2(3, it)] = t1;

		//double total_read = double(this->aligned_block.size()) * ops->it_batch * 2;// forward-backward and we read all
		double total_read = double(this->total_read) * ops->it_batch * 2;// forward-backward writes
		double total_write = double(this->total_write) * ops->it_batch * 2;// forward-backward writes
		double GBps = (total_read + total_write) / (t1 * 1024 * 1024 * 1024);

		debug::stream << "it" << total_it << " LB=" << txt::String::Format("%6.4f", LB) << " E=" << txt::String::Format("%6.4f", best_E) << " dPhi=" << txt::String::Format("%1.0g", dPhi) << txt::String::Format(" mrg=%4.4g", margin) << " dt=" << txt::String::Format("%4.3f", t1) << "s. "
		<< txt::String::Format("%3.1f", GBps) << " GB/s.\n";

		//reeturn with stopping conditions in this order: zero_gap, target_energy, precision
		if (best_E - LB < ops->gap_tol){
			debug::stream << "Reached duality gap < " << best_E - LB << "\n";
			exit_reason = zero_gap;
			//converged = true;
			break;
		};

		if (current_E < target_E){
			debug::stream << "Reached energy below " << target_E << "\n";
			exit_reason = target_energy;
			break;
		};

		if (dPhi < ops->conv_tol){
			debug::stream << "Converged to message precision " << dPhi << ".\n";
			//converged = true;
			exit_reason = precision;
			break;
		};

		if (LB < LB0){
			debug::stream << "Converged to LB precision " << (LB0 - LB)/LB0 << ".\n";
			//converged = true;
			exit_reason = precision;
			break;
		};
		LB0 = LB;
	};
	//total_it += it;
	if (exit_reason == error){
		debug::stream << "Reached maximum number of iterations (" << maxit << ").\n";
		exit_reason = iterations;
	};
	if (it < Max_it)++it; // number of iterations actually performed 
	debug::stream << "LB= " << txt::String::Format("%8.4f", LB) << " E= " << txt::String::Format("%8.4f", best_E) << "\n";
	//debug::stream << "TRWS Time: " << txt::String::Format("%3.2f", total_t1) << " s. ";
	hist.resize(mint2(4, it)); // in case when finished earlier
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::destroy_core(){
	aligned_block.destroy();
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::init(energy<type> * E){
	debug::PerformanceCounter c1;
	this->E = E;
	ops->conv_tol = E->tolerance*ops->rel_conv_tol;
	ops->gap_tol = E->tolerance*ops->rel_gap_tol;
	SG.init(E->G);
	debug::stream << "Stream graph layers: " << SG.nn.size() << " time: " << c1.time() << "s. \n";
	this->nV = E->nV();
	this->nE = E->nE();
	this->maxK = E->maxK;
	current_x.resize(nV);
	current_x << 0;
	best_x.resize(nV);
	best_x << 0;
	nodes.resize(nV);
	edges.resize(nE);
	init_core();
};

template<class type, class vectorizer>
term2v<type, vectorizer> * trws_machine<type, vectorizer>::construct_f2(int e, aallocator * al){
	// default to construct from E0
	return dynamic_cast<t_f2 *>(E->f2(e).copy(al));
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::init_core(){
	total_write = 0;
	// calculate size and placement, then allocate
	for (full_iterator ii(this); ii.valid(); ++ii){
		int s = ii.s();
		//debug::stream << s << " : " << aligned_block.mem_top - aligned_block.mem_beg << "\n";
		//};
		//for (int l = 0; l < SG.nn.size(); ++l){
		//for (int i = 0; i < SG.nn[l]; ++i){
		//int s = SG.levels[SG.level_start[l] + i];
		int K = E->K(s);
		// node.core 
		aligned_block.allocate_a<node_info_core>();
		// edge_info_core in_list
		aligned_block.allocate<edge_info_core>(E->G.in[s].size());
		// edge_info_core out_list
		aligned_block.allocate<edge_info_core>(E->G.out[s].size());
		// data
		int V = sizeof(vectorizer) / sizeof(type);
		int vK = (K + V - 1) / V; // round up to to a full vector
		// calculate required block size: f1 + bw + nout*msg
		aligned_block.allocate_a<vectorizer>(vK * 2);
		//
		//aligned_block.allocate_a<vectorizer>(vK);
		//
		total_write += (vK * sizeof(vectorizer)); // writes bw
		for (int j = 0; j < E->G.out[s].size(); ++j){
			int e = E->G.out[s][j];
			int t = E->G.E[e][1];
			int K2 = E->K(t);
			int vK2 = (K2 + V - 1) / V; // round up to to a full vector
			int mvK = std::max(vK, vK2);
			aligned_block.allocate_a<vectorizer>(mvK); // messages, compacted
			total_write += (vK + vK2) * sizeof(vectorizer) / 2; // average written messages
			// extra data per edge
			//edges[e].f2()->copy(&aligned_block);
			//aligned_block.allocate_a<vectorizer>(vK);
			//aligned_block.allocate_a<vectorizer>(vK2);
		};
		// extra data for f2: per edge, after messages
		for (int j = 0; j < E->G.out[s].size(); ++j){
			int e = E->G.out[s][j];
			construct_f2(e, &aligned_block);
		};
		//		};
	};
	// allocate r1 in the end (essentially out of core)
	for (full_iterator ii(this); ii.valid(); ++ii){
		int s = ii.s();
		int K = E->K(s);
		int V = sizeof(vectorizer) / sizeof(type);
		int vK = (K + V - 1) / V;
		aligned_block.allocate_a<vectorizer>(vK);
	};
	aligned_block.init();
	debug::stream << "core size:" << aligned_block.size() / 1024 / 1024 << "Mb, writable: " << total_write / 1024 / 1024 << "Mb.\n";
	// nodes
	//for (int l = 0; l < SG.nn.size(); ++l){
	//for (int i = 0; i < SG.nn[l]; ++i){
	//int s = SG.levels[SG.level_start[l] + i];
	for (full_iterator ii(this); ii.valid(); ++ii){
		int s = ii.s();
		node_info & node = nodes[s];
		int K = E->K(s);
		// node_info->active
		node.active = false;
		// node_info->K
		node.K = K;
		// node_info->in
		node.in.resize(E->G.in[s].size());
		// node_info->out
		node.out.resize(E->G.out[s].size());
		//node_info_core
		node.core = aligned_block.allocate_a<node_info_core>();
		node.core->n_in = E->G.in[s].size();
		node.core->n_out = E->G.out[s].size();
		// edge_info_core in_list, the pointer is calculated from node.core
		aligned_block.allocate<edge_info_core>(E->G.in[s].size());
		// edge_info_core out_list
		aligned_block.allocate<edge_info_core>(E->G.out[s].size());
		//
		// node_core_info->vK
		int V = sizeof(vectorizer) / sizeof(type);
		int vK = (K + V - 1) / V; // round up to to a full vector
		node.core->vK = vK;
		// node_info->core->node
		node.core->node = &node;
		//node.core->data: bw + f1
		node.core->data = aligned_block.allocate_a<vectorizer>(vK * 2);
		
		// node_info->core->r1, node_info->r1
		// node_info->core->bw, node_info->bw
		node.core->_bw()[vK - 1] = 0; // padd with zeros
		node.bw.set_ref(node.core->_bw(), K);
		node.bw << 0;
		// node_info->core->f1, node_info->f1
		node.core->_f1()[vK - 1] = INF(type); // padd with infinity
		node.f1.set_ref(node.core->_f1(), K);
		node.f1 << E->f1[s];
		/*
		node.core->__r1 = aligned_block.allocate_a<vectorizer>(vK);
		node.r1.set_ref(node.core->r1(), K);
		node.r1 << 0;
		*/
		//
		// create edges
		for (int j = 0; j < node.out.size(); ++j){
			int e = E->G.out[s][j];
			node.out[j] = &edges[e];
			int t = E->G.E[e][1];
			int K2 = E->K(t);
			int vK2 = (K2 + V - 1) / V; // round up to to a full vector
			int vmK = std::max(vK, vK2);
			// edge_info->core
			edge_info & edge = edges[e];
			edge.core = &node.core->_out_begin()[j];
			// edge_info->core->msg
			edge.core->msg = aligned_block.allocate_a<vectorizer>(vmK);
			// edge_info->msg
			edge.core->msg[vmK - 1] = 0;
			edge.msg.set_ref(edge.core->msg, std::max(node.K, K2));
			edge.msg << 0;
			// edge_info->tail
			edge.tail = &node;
		};
		// edge_info->core->f2
		for (int j = 0; j < node.out.size(); ++j){
			int e = E->G.out[s][j];
			edge_info & edge = edges[e];
			t_f2 * f2 = construct_f2(e, &aligned_block);
			edge.core->f2 = f2;
		};

		for (int j = 0; j < node.in.size(); ++j){
			int e = E->G.in[s][j];
			// node_info->in[j]
			node.in[j] = &edges[e]; // initialized in out list
			// edge_info->head
			edges[e].head = &node;
		};
		// node_info->core->coeff
		node.core->coeff = vertex_coeff(node);
		node.core->init();
		int z = 0;
		//};
	};
	total_read = aligned_block.size();
	// allocate r1 in the end (essentially out of core)
	for (full_iterator ii(this); ii.valid(); ++ii){
		int s = ii.s();
		node_info & node = nodes[s];
		int K = E->K(s);
		int V = sizeof(vectorizer) / sizeof(type);
		int vK = (K + V - 1) / V;
		//		
		node.core->__r1 = aligned_block.allocate_a<vectorizer>(vK);
		node.r1.set_ref(node.core->r1(), K);
		node.r1 << 0;
	};

	// Thi linkage pass
	// pass once again, fill: next, in_msg, in_f2
	//for (int l = 0; l < SG.nn.size(); ++l){ //parallel inside layer
	//for (int i = 0; i < SG.nn[l]; ++i){
	//int s = SG.levels[SG.level_start[l] + i];
	for (full_iterator ii(this); ii.valid(); ++ii){
		int s = ii.s();
		node_info & node = nodes[s];
		// core_node_info->next
		//if (i < SG.nn[l] - 1){
		//int t = SG.levels[SG.level_start[l] + i + 1];
		if (ii.next().valid()){
			int t = ii.next().s();
			node.core->next = nodes[t].core;
		} else{
			node.core->next = 0;
		};

		// incoming
		for (int j = 0; j < node.in.size(); ++j){
			edge_info * edge = node.in[j];
			edge_info_core * e_core = &node.core->_in_begin()[j];
			// edge_info->core->f2
			e_core->f2 = edge->core->f2;
			// edge_info->core->msg
			e_core->msg = edge->core->msg;
			// edge_info->core->_head
			e_core->_head = edge->head->core;
			e_core->_tail = edge->tail->core;
		};

		// outcoming
		for (int j = 0; j < node.out.size(); ++j){
			edge_info * edge = node.out[j];
			edge_info_core * e_core = &node.core->out_begin()[j];
			// edge_info->core->_head, tail
			e_core->_head = edge->head->core;
			e_core->_tail = edge->tail->core;
		};
	};
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::set_f2(edge_info & ee, t_f2 * f2){
	ee.core->f2 = f2;
	//ee.core->f2_data = f2->data_ptr();
	node_info * tail = ee.tail;
	int j1 = tail->find_out(&ee);
	{
		edge_info_core & ec = tail->core->out_begin()[j1];
		ec.f2 = f2;
		//ec.f2_data = f2->data_ptr();
	};
	node_info * head = ee.head;
	int j2 = head->find_in(&ee);
	{
		edge_info_core & ec = head->core->in_begin()[j2];
		ec.f2 = f2;
		//ec.f2_data = f2->data_ptr();
	};
};

template<class type, class vectorizer>
double trws_machine<type, vectorizer>::cost(const intf & x)const{
	double r = 0;
	for (int s = 0; s < nV; ++s){
		const node_info & v = nodes[s];
		r += v.f1[x[s]];
		for (int j = 0; j < v.out.size(); ++j){
			const edge_info & edge = *v.out[j];
			int t = int(edge.head - nodes.begin());
			int x_s = x[s];
			int x_t = x[t];
			if(true){//if (!edge.transposed){
				r += edge.core->f2->operator()(x_s, x_t);
			} else{
				r += edge.core->f2->operator()(x_t, x_s);
			};
		};
	};
	//debug::stream << "E.cost = " << E->cost(x) << "\n";
	return r;
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::set_M(const num_array<double, 2> & M){
	if (M.size()[1] < nE){
		debug::stream << "initial M is empty\n";
		for (int e = 0; e < nE; ++e){
			edges[e].msg << 0;
		};
		return; // empty input messages
	};
	for (int e = 0; e < nE; ++e){
		edges[e].msg << M.subdim<1>(e);
	};
	init_iteration(); // restore bw sums
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::get_M(num_array<double, 2> & M){
	//M.resize(mint2(K, nE));
	for (int e = 0; e < nE; ++e){
		M.subdim<1>(e) << edges[e].msg;
		/*
		for (int i = 0; i < M.size()[0]; ++i){
			M(e, i) = edges[e].msg[i];
		};
		*/
	};
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::get_message(int e, bool forward, std::vector<type> & x){
	if (!forward){
		edge_info & ee = edges[e];
		node_info * S = ee.tail;
		node_info * T = ee.head;
		x.resize(S->K);
		for (int i = 0; i < S->K; ++i){
			x[i] = ee.msg[i];
		};
	} else{
		edge_info & ee = edges[e];
		node_info * S = ee.tail;
		node_info * T = ee.head;
		x.resize(T->K);
		for (int i = 0; i < T->K; ++i){
			type v = 0;
			if (T->f1[i] < INF(type)){
				v = (T->r1[i] / T->core->coeff - T->f1[i] - T->bw[i]) / T->in.size();
			};
			x[i] = v;
		};
		/*
		tvect m(x.begin(), K);
		tvect temp = S->r1;
		temp -= ee.msg;
		ee.f2->min_sum(temp,m);
		*/
	};
};

/*
template<class type, class vectorizer> void trws_machine<type, vectorizer>::f_msgs::init(){
	msg_valid = false;
};

template<class type, class vectorizer> void trws_machine<type, vectorizer>::f_msgs::update(){
	m.resize(p->nE);
	// forward pass without averaging
	for (full_iterator ii(this); ii.valid(); ++ii){
		int s = ii.s();
		node_info & v = nodes[s];
		for (int j = 0; j < v.out.size(); ++j){
				edge_info * e = v.out[j];
			};
		};
	};
};

template<class type, class vectorizer> void trws_machine<type, vectorizer>::f_msgs::invalidate(){
	msg_valid = false;
};
*/

template<class type, class vectorizer>
type trws_machine<type, vectorizer>::vertex_coeff(node_info & v, type multiplier){
	type coeff;
	int n_in = v.in.size();
	int n_out = v.out.size();
	if (n_in == 0 && n_out == 0){
		return 1;
	};
	if (multiplier > 0){
		coeff = type(multiplier*1.0 / std::max(n_in, n_out));
	} else{
		coeff = -multiplier;
	};
	return coeff;
};

template<class type, class vectorizer>
trws_machine<type, vectorizer>::trws_machine(){
	//info_core = 0;
	ops = &_ops;
//	debug_LB = false;
};

template<class type, class vectorizer>
trws_machine<type, vectorizer>::~trws_machine(){
	//if (info_core){
	//	mfree(info_core);
	//};
};

template<class type, class vectorizer>
void trws_machine<type, vectorizer>::check_get_col(){
#ifdef _DEBUG
	const int V = sizeof(vectorizer) / sizeof(type);
	int maxKa = v_align<type, vectorizer>(maxK);
	tvect t2(maxKa);
	dynamic::fixed_vect<type, array_allocator<type, 16> > t1(maxKa);
	vectorizer * v = (vectorizer*)t1.begin();
	for (int e = 0; e < nE; ++e){
		int s = E->G.E[e][0];
		int t = E->G.E[e][1];
		node_info & T = nodes[t];
		int y_s = nodes[s].core->y;
		int y_t = nodes[t].core->y;
		int Ka = v_align<type,vectorizer>(nodes[s].K);
		tvect t1r;
		t1r.set_ref(v, Ka);
		for (int j0 = 0; j0 < nodes[t].K; ++j0){
			if (T.f1.begin()[j0] < INF(type)){ //get_col is only for movable and y
				t2.resize(Ka);
				edges[e].f2()->get_col(j0, t2);
				edges[e].f2()->get_col_e(nodes[s].core->vK, j0, v);
				for (int i = 0; i < Ka; ++i){
					node_info & S = nodes[s];
					if (S.f1.begin()[i] < INF(type)){
						type v1 = t2[i];
						type v2 = t1r[i];
						if (std::abs(v1 - v2) > 1e-4){
							edges[e].f2()->get_col_e(nodes[s].core->vK, j0, v);
						};
					};
				};
			};
		};
	};
#endif
};
template<typename type, typename vectorizer>
int trws_machine<type, vectorizer>::get_num_chunks()const{
	int NT = 1;
#if defined(USE_OPENMP) // || defined(_OPENMP)	
	NT = ops->max_CPU;
#endif
	return NT;
};

template<typename type, typename vectorizer>
bool trws_machine<type, vectorizer>::chunk_iterator::valid()const{
	return th_id < p->get_num_chunks();
};

template<typename type, typename vectorizer>
bool trws_machine<type, vectorizer>::slice_iterator::valid()const{
	if (!ci.valid())return false;
	return(l < ci.p->SG.nn.size());
};

template<typename type, typename vectorizer>
int trws_machine<type, vectorizer>::slice_iterator::line_begin()const{
	if (!valid())return -1;
	if (th_id() == 0){
		return 0;
	};
	/*
	int n = ci.p->SG.nn[l]; // number of nodes in line
	int NT = ci.p->get_num_chunks();
	int block = n / NT;
	return block*th_id() + (n%NT);
	*/
	int n = ci.p->SG.nn[l]; // number of nodes in line
	int NT = ci.p->get_num_chunks(); // number chunks
	int block = n / NT;
	if (n%NT > 0)block += 1;
	return block*th_id();
};

template<typename type, typename vectorizer>
int trws_machine<type, vectorizer>::slice_iterator::line_end()const{
	if (!valid())return -1;
	/*
	int n = ci.p->SG.nn[l]; // number of nodes in line
	int NT = ci.p->get_num_chunks();
	int block = n / NT;
	return block*(th_id()+1) + (n%NT);
	*/
	int n = ci.p->SG.nn[l]; // number of nodes in line
	int NT = ci.p->get_num_chunks(); // number chunks
	int block = n / NT;
	if (n%NT > 0)block += 1;
	int t = th_id();
	int end = block*(t+1);
	if (end > n){
		end = n;
	};
	return end;
};

template<typename type, typename vectorizer>
typename trws_machine<type, vectorizer>::node_info * trws_machine<type, vectorizer>::slice_iterator::first()const{
	line_iterator li(*this);
	if (!li.valid()){
		return 0;
	};
	return &ci.p->nodes[li.s()];
};

template<typename type, typename vectorizer>
typename trws_machine<type, vectorizer>::node_info * trws_machine<type, vectorizer>::slice_iterator::last()const{
	int i = line_end() - 1;
	if (i < line_begin()){
		return 0;
	};
	int s = ci.p->SG.levels[ci.p->SG.level_start[l] + i];
	return &ci.p->nodes[s];
};

template<typename type, typename vectorizer>
bool trws_machine<type, vectorizer>::line_iterator::valid()const{
	if (!si.valid())return false;
	return i < si.line_end();
};

template<typename type, typename vectorizer>
int trws_machine<type, vectorizer>::line_iterator::s()const{
	return p->SG.levels[p->SG.level_start[l()] + i];
};

//_____________instances______________________
template class trws_machine < float, float >;
template class trws_machine < double, double >;
template class trws_machine < float, sse_float_4 >;