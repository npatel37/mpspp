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
typedef PsimagLite::Geometry<RealType,InputNgType::Readable,Mpspp::ProgramGlobals> GeometryType;
typedef Mpspp::ParametersMpsSolver<RealType,ComplexOrRealType,InputNgType::Readable>
ParametersSolverType;
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
			  const ModelBaseType& model,
			  InputValidatorType& io)
{
	Mpspp::MpsSolver<ModelBaseType,InternalProductTemplate> mpsSolver(mpsSolverParams,
																	  model,
																	  io);
	mpsSolver.computeGroundState();
	std::cout<<"That's all folks!\n";
}

int main(int argc,char *argv[])
{
	Mpspp::InputCheck inputCheck;
	PsimagLite::String filename="";
	int opt = 0;
	PsimagLite::String strUsage(argv[0]);
	strUsage += " -f filename";
	int precision = 6;
	bool versionOnly = false;

	while ((opt = getopt(argc, argv,"f:p:V")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 'p':
			precision = atoi(optarg);
			std::cout.precision(precision);
			std::cerr.precision(precision);
			break;
		case 'V':
			versionOnly = true;
			break;
		default:
			inputCheck.usageMain(strUsage);
			return 1;
		}
	}

	// sanity checks here
	if (filename=="" && !versionOnly) {
		inputCheck.usageMain(strUsage);
		return 1;
	}

	//! setup distributed parallelization
	SizeType npthreads = 1;
	ConcurrencyType concurrency(&argc,&argv,npthreads);

	// print license
	if (ConcurrencyType::root()) {
		std::cout<<Mpspp::ProgramGlobals::license;
		Provenance provenance;
		std::cout<<provenance;
	}

	if (versionOnly) return 0;

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryType geometry(io);

	ParametersSolverType mpsSolverParams(io);

	ModelSelectorType modelSelector(mpsSolverParams.model);

	const ModelBaseType& model = modelSelector(mpsSolverParams,io,geometry);

	if (mpsSolverParams.options.find("InternalProductStored")!=PsimagLite::String::npos) {
		mainLoop<ModelBaseType,Mpspp::InternalProductStored>(mpsSolverParams,model,io);
	} else {
		mainLoop<ModelBaseType,Mpspp::InternalProductOnTheFly>(mpsSolverParams,model,io);
	}
}

/*@}*/

