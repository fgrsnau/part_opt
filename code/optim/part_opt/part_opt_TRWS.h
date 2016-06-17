#ifndef part_opt_trws_h
#define part_opt_trws_h

#include <queue>


#include "dynamic/num_array.h"
#include "dynamic/array_allocator.h"
#include "dynamic/fixed_array1.h"
#include "dynamic/fixed_array2.h"
#include "exttype/fixed_vect.h"
#include "optim/graph/mgraph.h"
#include "optim/trws/stream_graph.h"
#include "debug/performance.h"

#include "energy.h"
#include "trws_machine.h"

using namespace dynamic;

#define FINF std::numeric_limits<VALUE_TYPE>::infinity();
//typedef sse_float_4 vectorizer;
//typedef d_type type;
//typedef d_vectorizer vectorizer;
//typedef float_v4 vtype;
//typedef trws_machine<vtype> alg_trws;
//typedef energy<float> tenergy;
//
//_________________________alg_po_trws_________________________________
//! class for maximum persistency with message passing solver
template<class vtype = float_v4>
class alg_po_trws : public trws_machine<vtype>{
public:
	typedef trws_machine<vtype> parent;
	typedef typename trws_machine<vtype>::t_f2 t_f2;
	typedef energy<typename vtype::type> tenergy;
	typedef typename parent::node_info node_info;
	typedef typename parent::edge_info edge_info;
	typedef typename parent::tvect tvect;
	typedef typename parent::type type;
	using parent::nV;
	using parent::nE;
	using parent::maxK;
	using parent::nodes;
	using parent::edges;
	using parent::best_E;
	using parent::best_x;
	using parent::total_it;
	using parent::cost;
	using parent::destroy_core;
	using parent::init_core;
	using parent::init_iteration;
	using parent::exit_reason;
	using parent::check_get_col;
	using parent::run_converge;
	using parent::zero_gap;
	using parent::target_energy;
public:
	tenergy * E0; // initial input energy
	//energy F; // manipulated energy
	//dynamic::fixed_array1<t_f2*> ff2;
public: //data
	fixed_array1<po_mask> UU; // lists of immovable nodes
	std::queue<node_info*> dee_q;
public: //exported data
	num_array<int, 2> P; // improving mapping
	num_array<double, 3> * burn;
public: // options
	alg_po_trws_ops ops;
	intf y; // labeling to which we will reduce
private://internal and temp vars
	tvect v1;
	tvect v2;
	bool f2_replaced;
	intf y1;
private:
	static intf y0;
	int dee_mark;
	int dee_mark_strong;
public://statistics
	int nimmovable;
	//int total_it;
	int maxelim;
	int po_it;
private:
	virtual t_f2 * construct_f2(int e, aallocator * al)override;
protected:
	void finish();
	//void cut(intf & x);
	bool cut_tr(char * caller);
	bool flat_cuts();
	bool apply_cut(intf & x, char * caller);
	void set_y(intf & y);
	void rebuild_energy();
	void reduce_immovable();
	void rebuild_incremental(node_info & v, int k);
	//void expansion_move();
	double dee_delta(node_info & v, int k);
	double icm_delta(node_info & v, int k);
	void icm();
	int mark_dee();
	void mark_WTA();
	void unmark();
	void queue_neib(node_info & v);
	typedef enum { INIT = 0, WTA = 1, CUT = 2, PXCUT = 3} prunner;
	void mark_immovable_label(node_info & v, int k, prunner who);
	int mark_immovable(node_info & v, int k, prunner who);
	void report_po_iter(int total_mark1, int total_mark2, std::string caller);
	void init_aux();
	bool check_zero_gap();
	//
	std::pair<double,int> marg(int s);
	bool remove_condition(int s);
	bool stop_condition();
	void cleanup();
//==========================PUBLIC===============================================
public: //input
	alg_po_trws();
	~alg_po_trws();
	void set_E(tenergy * E);
	void init();
	void init1();
	void set_P(num_array<int, 2> & P); //!< export improving mapping
public: //run
	void prove_optimality(intf & y = y0);//!< Start PO iterations with target labeling y if given or find it by TRWS
public: //results:
	bool is_alive(int s, int k); //!< True if label k in node s is alive (not eliminated)
	int test_label(int s); //!< Returns the test labeling y in node s
public: //statistics
	double elim_labels_percent();
public: //performance
	debug::PerformanceCounter po_cuts;
	debug::PerformanceCounter po_msgs;
	debug::PerformanceCounter po_node_cuts;
	debug::PerformanceCounter po_total;
	debug::PerformanceCounter po_rebuild;
//===============================================================================
};

#endif