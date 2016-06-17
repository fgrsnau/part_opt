#ifndef energy_h
#define energy_h

#include "dynamic/fixed_array1.h"
#include "dynamic/fixed_array2.h"
#include "exttype/fixed_vect.h"
#include "optim/graph/mgraph.h"
#include "geom/math.h"
#include <bitset>
#include <limits>

#include "vectorizers.h"
#include "aallocator.h"

#include <typeinfo>
#include <typeindex>

#define INF(type) std::numeric_limits<type>::infinity()
typedef std::bitset<512> po_mask;

using namespace exttype;
using namespace dynamic;

// aligns size measured in types to integer # of vectorizers
template<typename type, typename vectorizer> int v_align(int size);
template<class vtype> int v_align(int size);
template<class vtype> mint2 v_align(mint2 size);


// itroduce energy with virtual pairwise functions, currently a virtual function returns full pw term
// abstract pw function:
// needed capablities: 
// min_sum, min_sum_t
// get_row, get_column, get_full
// min_sum_mask, min_sum_mask_t

class term2_base{
public:
	virtual ~term2_base(){};
	virtual term2_base * construct_copy(aallocator * al, std::type_index id) const = 0;
};

template<typename type>
class term2 : public term2_base{
public:
	typedef exttype::fixed_vect<type> tvect;
	typedef term2<type> tthis;
public:
	virtual ~term2(){};
	virtual void min_sum(tvect & a, tvect & r)const = 0;
	virtual void min_sum_t(tvect & a, tvect & r)const{
		return this->min_sum(a, r);
	};
	virtual type operator()(int i, int j) const = 0;
	virtual void get_row(int i, tvect & r) const = 0;
	virtual void get_col(int j, tvect & r) const = 0;
public: // some compatibility stuff
	virtual int count(int i) const = 0;
public: // methods for partial optimality
	virtual void reduce(po_mask & U_s, int y_s, po_mask & U_t, int y_t){};//defaults to empty
	virtual void rebuild(po_mask & U_s, int y_s, po_mask & U_t, int y_t){};//defaults to empty
public:
	//virtual tthis * copy(aallocator * al)const = 0;
	//virtual void * data_ptr(){return 0; };
};

template<typename vtype>
class term2v : public term2<typename vtype::type>{
public:
	typedef typename vtype::type type;
	typedef typename vtype::vectorizer vectorizer;
	typedef typename term2<type>::tvect tvect;
public: // optimized interface:
	virtual void get_col_e(const int vK1, int j0, vectorizer * message)const = 0;
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const{
		// this is to check correctness, otherwise nested virtual function is not good
		// default to min_sum regular vector implementation
		int V = sizeof(vectorizer) / sizeof(type);
		int K1 = vK1*V;
		tvect a; a.set_ref(source, K1);
		tvect r; r.set_ref(message, K1);
		this->min_sum(a, r);
	};
	virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const{
		// default to min_sum_t regular vector implementation
		int V = sizeof(vectorizer) / sizeof(type);
		int K1 = vK1*V;
		tvect a; a.set_ref(source, K1);
		tvect r; r.set_ref(message, K1);
		this->min_sum_t(a, r);
	};
};

template <template<class> class target, typename vtype> class term2v_CRTP_base : public term2v<vtype>{
public:
	typedef target<vtype> tsrc;
	tsrc & self(){ return dynamic_cast<tsrc &>(*this); };
	const tsrc & self() const { return dynamic_cast<const tsrc &>(*this); };
public:
	virtual term2_base * construct_copy(aallocator * al, std::type_index id) const override{
		if (id == std::type_index(typeid(float_v1))){
			return al->allocate< target<float_v1>, tsrc >(self());
		};
		if (id == std::type_index(typeid(float_v4))){
			return al->allocate< target<float_v4>, tsrc >(self());
		};
		if (id == std::type_index(typeid(double_v1))){
			return al->allocate< target<double_v1>, tsrc >(self());
		};
		throw debug_exception("not a recognizable vectorizer");
	};
};
/*
// float
template <template<class, class> class target, typename vectorizer> class term2v_CRTP_base<target, float, vectorizer> : public term2v<float, vectorizer>{
public:
	typedef target<float, vectorizer> tsrc;
	tsrc & self(){ return dynamic_cast<tsrc &>(*this); };
	const tsrc & self() const { return dynamic_cast<const tsrc &>(*this); };
public:
	virtual term2<float> * construct_copy(aallocator * al, std::type_index id) const override{
		if (id == std::type_index(typeid(float))){
			return al->allocate< target<float, float>, tsrc >(self());
		};
		if (id == std::type_index(typeid(sse_float_4))){
			return al->allocate< target<float, sse_float_4>, tsrc >(self());
		};
		throw debug_exception("not a recognizable vectorizer");
	};
};

// double
template <template<class, class> class target, typename vectorizer> class term2v_CRTP_base<target, double, vectorizer> : public term2v<double, vectorizer>{
public:
	typedef target<double, vectorizer> tsrc;
	tsrc & self(){ return dynamic_cast<tsrc &>(*this); };
	const tsrc & self() const { return dynamic_cast<const tsrc &>(*this); };
public:
	virtual term2<double> * construct_copy(aallocator * al, std::type_index id) const override{
		if (id == std::type_index(typeid(double))){
			return al->allocate< target<double, double>, tsrc >(self());
		};
		//if (id == std::type_index(typeid(sse_double_4))){
		//	return al->allocate< target<double, sse_double_4>, tsrc >(self());
		//};
		throw debug_exception("not a recognizable vectorizer");
	};
};
*/

/*
template<typename ttype, typename tvectorizer>
class term2_vcore{
public:
	typedef ttype type;
	typedef tvectorizer vectorizer;
public:
	term2<type> * src; // pointer to source term2
	//typedef typename term2<type>::tvect tvect;
public: // optimized interface:
	//! pass message in forward direction: minimize over first index
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const;
	//! pass message in backward direction 
	virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const;
	//! callback for allocating data and initialization
	virtual void init(aallocator & al, bool mem_valid);
};
*/

//_______________________term2v_matrix______________________________

template<typename vtype>
class term2v_matrix : public term2v_CRTP_base<term2v_matrix, vtype>, public num_array<typename vtype::type, 2, array_allocator<typename vtype::type, 16> > {//public exttype::ivector<fixed_matrix<type>, dynamic::fixed_array2<type> >{
public:
	//typedef typename vtype::type type;
	//typedef typename vtype::vectorizer vectorizer;
private:
	typedef num_array<type, 2, array_allocator<type, 16> > parent;
	typedef typename term2v<vtype>::tvect tvect;
	typedef term2v_matrix<vtype> tthis;
public:
	// virtual void * data_ptr()override{return this->begin(); };
public:
	using parent::begin;
	//empty constructor
	term2v_matrix(){};
	//copy with allocator
	template<class vtype2>
	term2v_matrix(const term2v_matrix<vtype2> & X, aallocator * al);
	//construct empty matrix of required size
	term2v_matrix(const mint2 & sz, aallocator * al);
	void construct(const mint2 & sz, aallocator * al);
	// construct from a num_array
	template<typename type2> explicit term2v_matrix(const num_array<type2, 2> & m);

	//virtual tthis * copy(aallocator * al)const override{
	//	return al->allocate<tthis>(*this);
	//};
	/*
	static tthis * allocate_copy(const term2<ttype> & x, aallocator * al)override{
		return al->allocate<tthis>(x);
	};
	*/
public: // implementing interfaces
	virtual type operator()(int i, int j) const override{
		return parent::operator()(i, j);
	};
	type & operator()(int i, int j){
		return parent::operator()(i, j);
	};
public:
	virtual void min_sum(tvect& src, tvect & dest)const override;
	virtual void min_sum_t(tvect& src, tvect& dest)const override;
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const override;
	virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const override;
	virtual void get_row(int i, tvect & r) const override{
		for (int j = 0; j < r.count(); ++j){
			r[j] = parent::operator()(i,j);
		};
	};
	virtual void get_col(int j, tvect & r) const override{
		assert(r.count() <= count(0));
		for (int i = 0; i < r.count(); ++i){
			r[i] = parent::operator()(i, j);
		};
	};
	virtual void get_col_e(const int vK1, int j0, vectorizer * message)const override{
		const vectorizer * pdata = (const vectorizer *)begin(0, j0);
		for (int i = 0; i < vK1; ++i){
			message[i] = pdata[i];
		};
	};
	virtual int count(int i) const override{
		return parent::size()[i];
	};
};

//_______________________term2v_diff______________________________

template<class vtype>
class term2v_diff : public term2v_CRTP_base<term2v_diff, vtype>{
private:
	typedef typename term2v<vtype>::tvect tvect;
	typedef term2v_diff<vtype> tthis;
//	template<class vtype2> friend class term2v_diff<vtype2>;
public:
	tvect diags; // assume matrix has constant diagonals
	tvect rdiags; // the same in the reverse order
	mint2 size; // aligned size
public:
	//virtual void * data_ptr()override{ return diags.begin(); };
public: // implementing interfaces
	virtual type operator()(int i, int j) const override{
		return diags(i - j + size[1] - 1);
	};
	type & operator()(int i, int j){
		return diags(i - j + size[1] - 1);
	};
public:
	//! construct from a different vectorizer
	template<class vtype2>
	term2v_diff(const term2v_diff<vtype2> & X, aallocator * al);
public:
	template<typename type2> explicit term2v_diff(const num_array<type2, 2> & m);
	virtual void min_sum(tvect& src, tvect & dest)const override;
	virtual void min_sum_t(tvect& src, tvect& dest)const override;
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const override;
	virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const override;
	virtual void get_row(int i, tvect & r) const override{
		for (int j = 0; j < size[1]; ++j){
			r[j] = term2v_diff::operator()(i, j);
		};
	};
	virtual void get_col(int j, tvect & r) const override{
		for (int i = 0; i < size[0]; ++i){
			r[i] = term2v_diff::operator()(i, j);
		};
	};
	virtual void get_col_e(const int vK1, int j0, vectorizer * message)const override;
	virtual int count(int i) const override{
		return size[i];
	};
};

//_______________________term2v_potts_______________________________
template<class vtype>
class term2v_potts : public term2v_CRTP_base<term2v_potts, vtype>{
	typedef typename term2v<vtype>::tvect tvect;
	typedef term2v_potts<vtype> tthis;
public:
	type gamma;
	mint2 size; // aligned size
public:
	template<typename type2> term2v_potts(const dynamic::num_array<type2, 2> & _f2);
	term2v_potts(mint2 sz, type gamma){
		if (gamma < 0)throw debug_exception("Will break down with negative gamma");
		this->gamma = gamma;
		this->size = sz;
	};
	term2v_potts(){};
	template<class vtype2>
	term2v_potts(const term2v_potts<vtype2> & a, aallocator * al = 0);
	//virtual tthis * copy(aallocator * al)const override;
	//
	virtual void min_sum(tvect & a, tvect & r)const override;
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const;
	virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const{
		term2v_potts::min_sum_e(vK1,source, message);
	};
	virtual type operator()(int i, int j)const override{
		if (i == j){
			//return type(-gamma);
			return 0;
		} else{
			//return 0;
			return gamma;
		};
	};
	virtual void get_row(int i, tvect & r) const override{
		//r << 0;
		//r[i] = -gamma;
		r << gamma;
		r[i] = 0;
	};
	virtual void get_col(int i, tvect & r) const override{
		r << gamma;
		r[i] = 0;
		//r << 0;
		//r[i] = -gamma;
	};
	virtual void get_col_e(const int vK1, int j0, vectorizer * message)const override{
		for (int i = 0; i < vK1; ++i){
			message[i] = gamma;
		};
		((type*)message)[j0] = 0;
	};
	virtual int count(int i) const override{
		return size[i];
	};
};

//_______________________term2v_tlinear_______________________________
template<class vtype>
class term2v_tlinear : public term2v_CRTP_base<term2v_tlinear, vtype>{
	// pairwise function min(gamma*|i-j|,th);
	typedef typename term2v<vtype>::tvect tvect;
	typedef term2v_tlinear<vtype> tthis;
public:
	type gamma;
	type th;
	mint2 size; // aligned size
public:
	term2v_tlinear(mint2 sz, type gamma, type th){
		if (gamma < 0)throw debug_exception("Will break down with negative gamma - non-convex");
		this->gamma = gamma;
		this->th = th;
		this->size = sz;
	};
	//! copy constructor
	template<class vtype2>
	term2v_tlinear(const term2v_tlinear<vtype2> & X, aallocator * al = 0){
		this->gamma = (type)X.gamma;
		this->th = (type)X.th;
		this->size = v_align<vtype>(X.size);
	};
	/*
	virtual tthis * copy(aallocator * al)const override{
		return al->allocate<tthis>(*this);
	};
	*/
	template<typename type2> term2v_tlinear(const dynamic::num_array<type2, 2> & _f2);
	virtual void min_sum(tvect & a, tvect & r)const override;
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const;
	__forceinline virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const{
		term2v_tlinear::min_sum_e(vK1, source, message);
	};
	virtual type operator()(int i, int j)const override{
		return std::min(gamma*abs(i-j),th);
	};
	virtual void get_row(int i, tvect & r) const override{
		for (int j = 0; j < size[1]; ++j){
			r[j] = term2v_tlinear::operator()(i, j);
		};
	};
	virtual void get_col(int j, tvect & r) const override{
		for (int i = 0; i < size[0]; ++i){
			r[i] = term2v_tlinear::operator()(i, j);
		};
	};
	virtual void get_col_e(const int vK1, int j0, vectorizer * message)const override{
		int i = 0;
		type * m = (type*)message;
		type v = gamma*j0;
		for (; i <= j0; ++i){
			m[i] = std::min(v, th);
			v -= gamma;
		};
		v  = gamma;
		for (; i < size[0]; ++i){
			m[i] = std::min(v, th);
			v += gamma;
		};
	};
	virtual int count(int i) const override{
		return size[i];
	};
};

//_______________________term2v_tquadratic_______________________________
template<class vtype>
class term2v_tquadratic : public term2v_CRTP_base<term2v_tquadratic, vtype>{
	// pairwise function min(gamma*|i-j|,th);
	typedef typename term2v<vtype>::tvect tvect;
	typedef term2v_tquadratic<vtype> tthis;
	struct q_roof{
		int i;  // when entered
		type a; // value a[i]
		int x; // when becomes active
	};
public:
	type gamma;
	type th;
	mint2 size; // aligned size
public:
	//! Try construct from a matrix and throw if not truncated quadratic
	template<typename type2> term2v_tquadratic(const num_array<type2,2> & matrix, mint2 sza);
	term2v_tquadratic(mint2 sz, type gamma, type th){
		if (gamma < 0)throw debug_exception("Will break down with negative gamma - non-convex");
		this->gamma = gamma;
		this->th = th;
		this->size = v_align<vtype>(sz);
	};
	template<class vtype2>
	term2v_tquadratic(const term2v_tquadratic<vtype2> & X, aallocator * al){
		this->gamma = (type)X.gamma;
		this->th = (type)X.th;
		this->size = v_align<vtype>(X.size);
	};
	/*
	virtual tthis * copy(aallocator * al)const override{
		return al->allocate<tthis>(*this);
	};
	*/
	virtual void min_sum(tvect & a, tvect & r)const override;
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const;
	__forceinline virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const{
		term2v_tquadratic::min_sum_e(vK1, source, message);
	};
	virtual type operator()(int i, int j)const override{
		return std::min(gamma*math::sqr(i - j), th);
	};
	virtual void get_row(int i, tvect & r) const override{
		for (int j = 0; j < size[1]; ++j){
			r[j] = term2v_tquadratic::operator()(i, j);
		};
	};
	virtual void get_col(int j, tvect & r) const override{
		for (int i = 0; i < size[0]; ++i){
			r[i] = term2v_tquadratic::operator()(i, j);
		};
	};
	virtual void get_col_e(const int vK1, int j0, vectorizer * message)const override{
		type *m = (type*)message;
		for (int i = 0; i < size[0]; ++i){
			m[i] = term2v_tquadratic::operator()(i, j0);
		};
	};
	virtual int count(int i) const override{
		return size[i];
	};
};


//_______________________term2v_matrix_po_____________________________
template<class vtype>
class term2v_matrix_po : public term2v_matrix<vtype>{
private:
	typedef term2v_matrix<vtype> parent;
public:
	const term2<type> * source;
public:
	term2v_matrix_po(){};
	~term2v_matrix_po(){
	};
	// construct from any term2v
	template<class vtype2>
	term2v_matrix_po(const term2v<vtype2> & m, aallocator * al) : parent(mint2(m.count(0), m.count(1)), al){
		source = &m;
	};
	virtual void reduce(po_mask & U_s, int y_s, po_mask & U_t, int y_t) override;
	virtual void rebuild(po_mask & U_s, int y_s, po_mask & U_t, int y_t) override;
};


//_______________________term2v_po_rduced_____________________
template<class vtype, template <typename> class base>
class term2v_po_reduced : public term2v<vtype>{
public:
	typedef typename vtype::type type;
	typedef typename vtype::vectorizer vectorizer;
	typedef base<vtype> tbase;
	typedef typename term2v<vtype>::tvect tvect;
	typedef term2v_po_reduced<vtype, base> tthis;
	tbase _src;
	tbase & src(){ return _src; };
	const tbase & src()const{ return _src; };
private:
	type c;
	vectorizer * _delta_st;
	vectorizer * _delta_ts;
	//dynamic::fixed_array1<vectorizer, array_allocator<vectorizer, 16> > block;
	int vK1;
	int vK2;
	//int vK1()const{ return base::count(0); };
	//int vK2()const{ return base::count(1); };
public:
	//virtual void * data_ptr()override{return _delta_st; };
	//virtual void * data_ptr()override{ return &src; };
public: //additional data
	tvect delta_st;
	tvect delta_ts;
	po_mask * U_s;
	po_mask * U_t;
#ifdef _DEBUG
	int y_s;
	int y_t;
	term2v_matrix_po<vtype> ref;
#endif
public:
	const tthis & self()const{ return *this; };
	virtual term2_base * construct_copy(aallocator * al, std::type_index id) const override{
		if (id == std::type_index(typeid(float_v1))){
			return al->allocate< term2v_po_reduced<float_v1, base> >(self());
		};
		if (id == std::type_index(typeid(float_v4))){
			return al->allocate< term2v_po_reduced<float_v4, base> >(self());
		};
		if (id == std::type_index(typeid(double_v1))){
			return al->allocate< term2v_po_reduced<double_v1, base> >(self());
		};
		throw debug_exception("not a recognizable vectorizer");
	};
	// copy constructor from base term2
	template<class vtype2>
	term2v_po_reduced(const base<vtype2> & Src, aallocator * al) :_src(Src, al){
		init(al);
	};
	// copy constructor from self
	template<class vtype2>
	term2v_po_reduced(const term2v_po_reduced<vtype2, base> & X, aallocator * al) :_src(X._src, al){
		init(al);
	};
	void init(aallocator * al);
	/*
	virtual tthis * copy(aallocator * al)const override{
		return al->allocate<tthis>(*this);
	};
	*/
public:
	virtual void reduce(po_mask & U_s, int y_s, po_mask & U_t, int y_t) override;
	virtual void rebuild(po_mask & U_s, int y_s, po_mask & U_t, int y_t) override;
public:
	virtual type operator()(int i, int j)const;
	type operator()(int i, int j, bool i_inU, bool j_inU)const;
	virtual void min_sum(tvect & a, tvect & r)const override;
	virtual void min_sum_t(tvect & a, tvect & r)const override;
public:
	template <bool dir_fw>
	void t_min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const;
	virtual void min_sum_e(const int vK1, vectorizer * source, vectorizer * message)const override;
	virtual void min_sum_et(const int vK1, vectorizer * source, vectorizer * message)const override;
private:
	void min_sum(tvect & a, tvect & r, po_mask * U_s, po_mask * U_t, const tvect & delta_st, const tvect & detla_ts)const;
	void check(tvect & a, tvect & r, tvect & r1, po_mask * U_s, po_mask * U_t, const tvect & delta_st, const tvect & detla_ts, int y_s, int y_t)const;
	virtual void get_row(int i, tvect & r)const override;
	virtual void get_col(int j, tvect & r)const override;
	virtual void get_col_e(const int vK1, int j0, vectorizer * message)const override;
	virtual int count(int i)const override{ return src().tbase::count(i); };
private: // instantiate some stuff
	void test(){
		base<float_v1> x1(_src, 0);
		base<float_v4> x2(_src, 0);
		base<double_v1> x3(_src, 0);
		base<double_v4> x4(_src, 0);
	};
};

//_______________________energy______________________________
template <class type>
class energy{
public:
	typedef typename scalalr_vectorizer<type>::vtype vtype;
	typedef typename term2<type>::tvect t_f1;
	typedef term2v<vtype> t_f2;
public:
	datastruct::mgraph G;
	intf K;
	int maxK;
	dynamic::fixed_array1<t_f1> f1;
	double tolerance;
public:
	virtual t_f2 & f2(int e) = 0;
	double cost(const intf & x);
	virtual void set_nE(int nE){
	};
	virtual void set_nV(int nV){
		G._nV = nV;
		f1.resize(nV);
		K.resize(nV);
	};
	
	void set_G(datastruct::mgraph & G){
		this->G = G;
		f1.resize(G.nV());
		K.resize(G.nV());
		set_nE(G.nE());
	};

	void set_E(dynamic::num_array<int, 2> E){
		int nE = E.size()[1];
		G.E.resize(nE);
		for (int e = 0; e < nE; ++ e){
			G.E[e][0] = E(0, e);
			G.E[e][1] = E(1, e);
		};
		set_nE(nE);
		G.edge_index();
	};
	int nV()const{
		return G.nV();
	};
	int nE()const{
		return G.nE();
	};
public:
	virtual ~energy(){};
public:
};


//____________________energy_auto_________________________________________

template<typename type>
class energy_auto : public energy<type>{
public:
	typedef typename energy<type>::t_f1 t_f1;
	typedef typename energy<type>::t_f2 t_f2;
public:
	dynamic::fixed_array1<t_f2*> F2;
public:
	int nfull;
	int npotts;
	int ntlinear;
	int ndiff;
	int ntquadratic;
	double maxf;
	double delta;
	double mult;
public:
	energy_auto();
	virtual void set_nE(int nE){
		this->G.E.resize(nE);
		this->F2.resize(nE);
	};
	virtual t_f2 & f2(int e){
		return *(this->F2[e]);
	};
	void cleanup();
	virtual ~energy_auto(){
		cleanup();
	};
	//bool test_potts(const dynamic::num_array<double, 2> & _f2);
	//bool test_tlinear(const dynamic::num_array<double, 2> & _f2);
	void set_f2(int e, const dynamic::num_array<double, 2> & _f2);
	void set_f1(int v, const dynamic::num_array<double, 1> & _f1);
public:
	void init();
	void report();
};


#endif