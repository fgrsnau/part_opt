#ifndef dee_h
#define dee_h

#include "dynamic/num_array.h"
#include "energy.h"

using namespace dynamic;
typedef double c_type;

double dee(const num_array<int,2>  & E, const num_array<c_type,2> & f1, const num_array<c_type,3> & f2, num_array<int,2> & X, num_array<int,2> & P,double eps=0);
double dee2(const num_array<int,2>  & E, const num_array<c_type,2> & f1, const num_array<c_type,3> & f2, num_array<int,2> & X, num_array<int,2> & P,double eps=0);

#endif