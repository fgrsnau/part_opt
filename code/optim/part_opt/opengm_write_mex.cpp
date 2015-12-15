#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/graphicalmodel/space/simplediscretespace.hxx>
#include <opengm/functions/potts.hxx>
#include <opengm/operations/adder.hxx>
#include <opengm/inference/messagepassing/messagepassing.hxx>
#include <opengm/inference/gibbs.hxx>
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#include <fstream>
#include <string>

#include "../src/interfaces/matlab/opengm/mex-src/model/matlabModelType.hxx"

#include "mex/mex_io.h"
#include "mex/mex_log.h"
#include "dynamic/num_array.h"
//#include "part_opt_TRWS.h"
#include "debug/performance.h"
#include "opengm_read.h"
#include "optim/part_opt/energy.h"


using namespace std; // 'using' is used only in example code
using namespace opengm;
using namespace exttype;


void mexFunction_protect(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	/*
	[] = opengm_write_mex(E,f1,f2,file_name);
	Input:
	E: [2 x nE] int32 - edges
	f1: [K x nV] double - unary costs
	f2: [K x K x nV] double - pairwise costs
	file_name string - output file
	X0: [K x nV] int32 - alive nodes
	*/
	if (nrhs < 4) {
		debug::stream << "Usage: [] = opengm_write_mex(E,f1,f2,file_name); \n";
		debug::stream << "Input:\n\
						 		\t E: [2 x nE] int32 - edges\n\
								\t f1 : [K x nV] double - unary costs\n\
								\t f2 : [K x K x nV] double - pairwise costs\n\
								\t file_name string - output file\n\
								\t X0 : [K x nV] int32 - alive mask\n";
		throw debug_exception("4 arguments expected\n");
	}

	mx_array<int, 2> E(prhs[0]);
	mx_array<double, 2> f1(prhs[1]);
	mx_array<double, 3> f2(prhs[2]);
	mx_string file_name(prhs[3]);
	
	int nV = f1.size()[1];
	int maxK = f1.size()[0];
	int nE = E.size()[1];

	num_array<int, 2> X0;
	if (nrhs == 5){
		X0 = mx_array<int, 2>(prhs[4]);
	} else{
		X0.resize(mint2(maxK, nV));
		X0 << 1;
	};

	//check sizes
	if (f2.size()[0] != maxK || f2.size()[1] != maxK || f2.size()[2] != nE){
		debug::stream << "f2 size expected " << maxK << "x" << maxK << "x" << nE << "\n";
		debug::stream << "provided" << f2.size() << "\n";
		throw debug_exception("wrong f2 size");
	};

	typedef double ValueType;
	/*typedef opengm::meta::TypeListGenerator
		<
		opengm::PottsFunction<ValueType>,
		opengm::ExplicitFunction<ValueType>
		> ::type FunctionTypeList;*/

	typedef opengm::meta::TypeListGenerator
		<
		opengm::PottsFunction<ValueType>,
		opengm::PottsNFunction<ValueType>,
		opengm::ExplicitFunction<ValueType>,
		opengm::TruncatedSquaredDifferenceFunction<ValueType>,
		opengm::TruncatedAbsoluteDifferenceFunction<ValueType>
		> ::type FunctionTypeList;

	typedef opengm::GraphicalModel<ValueType, opengm::Adder, FunctionTypeList>  GmType;

	typedef GmType::SpaceType Space;

	// count alive labels
	num_array<int, 1> K(nV);
	for (int u = 0; u < nV; ++u){
		int Ku = X0.subdim<1>(u).sum();
		K[u] = Ku;
	};
	Space space(K.begin(),K.end());

	GmType gm(space);

	//debug::stream << " # vars: " << gm.numberOfVariables() << "\n";
	debug::stream << "Model: nV = " << gm.numberOfVariables() << " maxK = " << maxK << " nE = " << nE << "\n";

	// construct graphical model
	// add unary terms

	for (int v = 0; v < nV; ++v){
		int Kv = K[v];
		const size_t shape[] = { Kv };
		ExplicitFunction<double> f(shape, shape + 1);
		int i1 = 0;
		for (int i = 0; i < maxK; ++i){
			if (!X0(i, v))continue;
			f(i1) = double(f1(i, v));
			++i1;
		};
		GmType::FunctionIdentifier fid = gm.addFunction(f);
		size_t variableIndices[] = { v };
		gm.addFactor(fid, variableIndices, variableIndices + 1);
	};
	// add pairwise terms
	int np = 0;
	int nf = 0;
	for (int e = 0; e < nE; ++e){
		int s = E(0, e);
		int t = E(1, e);
		size_t shape[] = { gm.numberOfLabels(s), gm.numberOfLabels(t) };

		GmType::FunctionIdentifier fid;
		// try potts
		try{
			term2v_potts<double> g(f2.subdim<2>(e));
			opengm::PottsFunction<double> f(gm.numberOfLabels(s), gm.numberOfLabels(t), 0, g.gamma);
			fid = gm.addFunction(f);
			++np;
		} catch (...){
			//then full
			opengm::ExplicitFunction<double> f(shape, shape + 2, 0.0);
			int i1 = 0;
			for (int k1 = 0; k1 < maxK; ++k1){
				if (!X0(k1, s))continue;
				int i2 = 0;
				for (int k2 = 0; k2 < maxK; ++k2){
					if (!X0(k2, t))continue;
					assert(i1 < K[s] && i2 < K[t]);
					f(i1, i2) = f2(k1, k2, e);
					++i2;
				};
				++i1;
			};
			fid = gm.addFunction(f);
			++nf;
		};

		size_t variableIndices[] = { s, t };
		//sort(variableIndices, variableIndices + 2);
		gm.addFactor(fid, variableIndices, variableIndices + 2);
	};
	debug::stream << "# Potts: " << np << "# Full: "<<nf<<"\n";
	debug::stream << "Saving " << file_name << "...";
	opengm::hdf5::save(gm, file_name, "gm");
	debug::stream << "ok.\n";
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	mexargs::MexLogStream log("log/output.txt", false);
	debug::stream.attach(&log);
	mexargs::MexLogStream err("log/errors.log", true);
	debug::errstream.attach(&err);

	//try{
	mexFunction_protect(nlhs, plhs, nrhs, prhs);
	//} catch (std::exception & e){
	//	debug::stream << "!Exception:" << e.what() << "\n";
	//};

	memserver::get_global()->clean_garbage();
	debug::stream.detach();
	debug::errstream.detach();
};