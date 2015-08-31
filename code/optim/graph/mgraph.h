#ifndef mgraph_h
#define mgraph_h

#include "dynamic/array_allocator.h"
#include "dynamic/fixed_array1.h"
#include "dynamic/fixed_array2.h"
#include "exttype/fixed_vect.h"
#include "dynamic/num_array.h"

namespace datastruct{
	using namespace exttype;
	using namespace dynamic;
	using exttype::mint2;

	class mgraph{
	public:
		typedef fixed_vect<int> intf;
		typedef dynamic::fixed_array1<mint2, dynamic::array_allocator<mint2> > tE;
		typedef dynamic::fixed_array1<intf, dynamic::array_allocator<intf > > tmy_vararray;
	public:
		//intf V;
		int _nV;
		tE E;
		tmy_vararray in;
		tmy_vararray out;
	public:
		mgraph(){
			_nV = 0;
		};
		void init(int nV,const dynamic::num_array<int,2> & E);
		//! number of nodes
		int nV()const{return _nV;};
		//int nV()const{return V.size();};
		//! number of edges
		int nE()const{return E.size();};
		//! number of edges adjucent to s (node degree)
		int nns(int s)const{return in[s].size()+out[s].size();};
	public:
		//! create a grid graph
		template<int rank, bool lastindexfastest> void create_grid(const intn<rank>& sz);
	public:
		//! initialize incoming and outcoming lists
		void edge_index();
	public:
		//! find edge s->t
		int edge(int s,int t);
		//! find edge t->s in the list of incoming edges to s
		int in_edge(int s,int t);
		//! find edge s->t in the list of outcoming edges from s
		int out_edge(int s,int t);
	public:
		//! Compute approximate minimul graph coloring
		int greedy_coloring(dynamic::fixed_array1<int> & coloring);
	public:
		//! Get edges as 2 x nE matrix of indices
		dynamic::num_array<int,2> get_E();
	};
};

#endif
