#ifndef vectorizers_h
#define vectorizers_h

#include "vectorizers1.h"

//class sse_double_4;

struct float_v1{
	typedef float type;
	typedef float vectorizer;
};

struct float_v4{
	typedef float type;
	typedef sse_float_4 vectorizer;
};

struct double_v1{
	typedef double type;
	typedef double vectorizer;
};

struct double_v4{
	typedef double type;
	typedef double vectorizer; // this is a stub, sse_double_4 is incomplete
};

template<typename type> struct  scalalr_vectorizer{
	typedef type vtype;
	typedef type vectorizer;
};
template<> struct  scalalr_vectorizer<float>{
	typedef float_v1 vtype;
	typedef float vectorizer;
};
template<> struct  scalalr_vectorizer<double>{
	typedef double_v1 vtype;
	typedef double vectorizer;
};

template<class type> class default_vectorizer{
public:
	typedef float_v1 vtype;
	typedef type vectorizer;
};

template<> class default_vectorizer<float>{
public:
	typedef float_v4 vtype;
	typedef float_v4::vectorizer vectorizer;
};

template<> class default_vectorizer<double>{
public:
	typedef double_v1 vtype;
	typedef double_v1::vectorizer vectorizer;
};


#endif