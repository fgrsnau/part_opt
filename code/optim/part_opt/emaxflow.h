#include "dynamic/num_array.h"
//
using namespace dynamic;

template<typename tcap> struct tflow{
	typedef tcap type;
};
template<> struct tflow<int>{
	typedef long long type;
};
template<> struct tflow<float>{
	typedef double type;
};

template<typename c_type>
double emaxflow(const num_array<int, 2>  & E, const num_array<c_type, 2> & f1, const num_array<c_type, 3> & f2, num_array<int, 1> & x);