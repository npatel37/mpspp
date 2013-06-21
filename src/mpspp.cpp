#include "String.h"
const PsimagLite::String license=
"Copyright (c) 2012, UT-Battelle, LLC\n"
"All rights reserved\n"
"[MPS++, Version 0.1]\n"
"\n"
"--------------------------------------------------------------------------------\n"
"\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED.\n"
"\n"
"Please see full disclaimer and open source license included in file LICENSE\n"
"--------------------------------------------------------------------------------\n"
"\n"
"\n";

/** \ingroup MPSPP */
/*@{*/

/*! \file mpspp.cpp
 *
 *  The MPS++ main driver
 *
 */

#include <unistd.h>
#include "ProgramGlobals.h"
#include "InputCheck.h"
#include "Provenance.h"
#include "Concurrency.h"
#include "InputNg.h"
#include "Geometry/Geometry.h"
#include "ParametersMpsSolver.h"
#include "ModelSelector.h"
#include "MpsSolver.h"
//#include "MpoLocal.h"
#include "InternalProductStored.h"
//#include "InternalProductKron.h"
#include "InternalProductOnTheFly.h"
#include "SymmetryLocal.h"

typedef double RealType;
typedef double ComplexOrRealType;
typedef PsimagLite::Concurrency ConcurrencyType;
typedef PsimagLite::InputNg<Mpspp::InputCheck> InputNgType;
typedef Mpspp::SymmetryLocal SymmetryLocalType;
typedef PsimagLite::Geometry<RealType,Mpspp::ProgramGlobals> GeometryType;
typedef Mpspp::ParametersMpsSolver<RealType,ComplexOrRealType,InputNgType::Readable> ParametersSolverType;
typedef PsimagLite::InputNg<Mpspp::InputCheck>::Readable InputValidatorType;

typedef Mpspp::ModelBase<ParametersSolverType,
						 InputValidatorType,
						 SymmetryLocalType,
						 GeometryType> ModelBaseType;
typedef Mpspp::ModelSelector<ModelBaseType> ModelSelectorType;
typedef ModelBaseType::MpoLocalType MpoLocalType;
typedef MpoLocalType::MpsLocalType MpsLocalType;

// FIXME: make configurable at runtime:

template<typename ModelBaseType,
template<typename,typename> class InternalProductTemplate>
void mainLoop(const typename ModelBaseType::ParametersSolverType& mpsSolverParams,
const ModelBaseType& model)
{
	Mpspp::MpsSolver<ModelBaseType,InternalProductTemplate> mpsSolver(mpsSolverParams,model);

	mpsSolver.computeGroundState();

//	const MpoLocalType& H = model.hamiltonian();

//	MpsLocalType hpsi = H*psi;

//	std::cout<<"Energy="<<scalarProduct(psi,hpsi)<<"\n";

	std::cout<<"That's all folks!\n";
}

int main(int argc,char *argv[])
{
	Mpspp::InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename";
	while ((opt = getopt(argc, argv,"f:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename=="") {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	//! setup distributed parallelization
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);

	// print license
	if (ConcurrencyType::root()) {
		std::cerr<<license;
		Provenance provenance;
		std::cout<<provenance;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryType geometry(io);

	ParametersSolverType mpsSolverParams(io);

	ModelSelectorType modelSelector(mpsSolverParams.model);

	const ModelBaseType& model = modelSelector(mpsSolverParams,io,geometry);

	//typename MpsLocalType::IoInputType ioForMps(mpsSolverParams.initialMps);

	//MpsLocalType psi(ioForMps);

	if (mpsSolverParams.options.find("InternalProductStored")!=PsimagLite::String::npos) {
		mainLoop<ModelBaseType,Mpspp::InternalProductStored>(mpsSolverParams,model);
		//} else if (mpsSolverParams.options.find("InternalProductKron")!=PsimagLite::String::npos) {
		//	mainLoop<ModelBaseType,Mpspp::InternalProductKron>(psi,mpsSolverParams,model);
	} else {
		mainLoop<ModelBaseType,Mpspp::InternalProductOnTheFly>(mpsSolverParams,model);
	}
}

/*@}*/

