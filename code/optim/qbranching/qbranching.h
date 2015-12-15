#ifndef q_branching_h
#define q_branching_h

#include "dynamic/options.h"
#include "dynamic/dynamic.h"

/*! qbranchin
  Input: ee: 2 x nE - Graph
         hh: 2 x nH - Line Graph
		 ww: 1 x nH - Line Graph costs
		 ops - optimization options
  Output: B: 2 x nV - solution
	      B(0,:) - incoming node index in ee, B(1,:) - incoming edge index in ee
*/

dynamic::num_array<int, 2> qbranching(dynamic::num_array<int, 2> & ee, dynamic::num_array<int, 2> & hh, dynamic::num_array<double, 1> & ww, options & ops);

#endif