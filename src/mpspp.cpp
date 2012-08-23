#include <string>
const std::string license=
"Copyright (c) 2012, UT-Battelle, LLC\n"
"All rights reserved\n"
"[MPS++, Version 0.1]\n"
"\n"
"********************************************************************************\n"
"DISCLAIMER\n"
"\n"
"THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND\n"
"CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED\n"
"WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n"
"PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE\n"
"COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,\n"
"OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR\n"
"ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR\n"
"CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,\n"
"PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
"DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\n"
"CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN\n"
"CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE\n"
"OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n"
"SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH\n"
"DAMAGE.\n"
"\n"
"NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED\n"
"STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR\n"
"ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY\n"
"INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS\n"
"DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.\n"
"\n"
"********************************************************************************\n"
"\n"
"Please see full open source license included in file LICENSE\n"
"\n";

/** \ingroup MPSPP */
/*@{*/

/*! \file mpspp.h
 *
 *  The MPS++ main driver
 *
 */

#include "ProgramGlobals.h"
#include "InputCheck.h"
#include "Provenance.h"
#include "ConcurrencySerial.h"
#include "InputNg.h"
#include "Geometry.h"
#include "ParametersMpsSolver.h"
#include "ModelSelector.h"
#include "MpsSolver.h"
#include "MatrixProductOperator.h"
#include "InternalProductStored.h"
//#include "InternalProductKron.h"
#include "InternalProductOnTheFly.h"
#include "SymmetryLocal.h"

typedef double RealType;
typedef double ComplexOrRealType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::InputNg<Mpspp::InputCheck> InputNgType;
typedef Mpspp::SymmetryLocal SymmetryLocalType;
typedef PsimagLite::Geometry<RealType,Mpspp::ProgramGlobals> GeometryType;
typedef Mpspp::ParametersMpsSolver<RealType,ComplexOrRealType,InputNgType::Readable> ParametersSolverType;
typedef PsimagLite::InputNg<Mpspp::InputCheck>::Readable InputValidatorType;

typedef Mpspp::ModelBase<ParametersSolverType,
						 InputValidatorType,
						 SymmetryLocalType,
						 GeometryType,
						 ConcurrencyType> ModelBaseType;
typedef Mpspp::ModelSelector<ModelBaseType> ModelSelectorType;
typedef ModelBaseType::MatrixProductOperatorType MatrixProductOperatorType;
typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;

// FIXME: make configurable at runtime:

template<typename ModelBaseType,
template<typename,typename> class InternalProductTemplate,
typename ConcurrencyType>
void mainLoop(typename ModelBaseType::MatrixProductOperatorType::MatrixProductStateType& psi,
const typename ModelBaseType::ParametersSolverType& mpsSolverParams,
const ModelBaseType& model,
ConcurrencyType& concurrency)
{
	Mpspp::MpsSolver<ModelBaseType,InternalProductTemplate> mpsSolver(mpsSolverParams,model,concurrency);

	mpsSolver.computeGroundState(psi);

//	const MatrixProductOperatorType& H = model.hamiltonian();

//	MatrixProductStateType hpsi = H*psi;

//	std::cout<<"Energy="<<scalarProduct(psi,hpsi)<<"\n";

//	std::cout<<"That's all folks!\n";
}

int main(int argc,char *argv[])
{
	Mpspp::InputCheck inputCheck;
	std::string filename="";
	int opt = 0;
	std::string strUsage(argv[0]);
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
	ConcurrencyType concurrency(argc,argv);

	// print license
	if (concurrency.root()) {
		std::cerr<<license;
		Provenance provenance;
		std::cout<<provenance;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	GeometryType geometry(io);

	ParametersSolverType mpsSolverParams(io);

	ModelSelectorType modelSelector(mpsSolverParams.model);

	const ModelBaseType& model = modelSelector(mpsSolverParams,io,geometry,concurrency);

	SymmetryLocalType symm;
	MatrixProductStateType psi(symm); // initialize to something

	if (mpsSolverParams.options.find("InternalProductOnTheFly")!=std::string::npos) {

		if (mpsSolverParams.options.find("InternalProductStored")!=std::string::npos) {
			mainLoop<ModelBaseType,Mpspp::InternalProductStored,ConcurrencyType>(psi,mpsSolverParams,model,concurrency);
			//} else if (mpsSolverParams.options.find("InternalProductKron")!=std::string::npos) {
			//	mainLoop<ModelBaseType,Mpspp::InternalProductKron,ConcurrencyType>(psi,mpsSolverParams,model,concurrency);
		} else {
			mainLoop<ModelBaseType,Mpspp::InternalProductOnTheFly,ConcurrencyType>(psi,mpsSolverParams,model,concurrency);
		}
	}
}

/*@}*/

