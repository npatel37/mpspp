#include <string>
const std::string license=
"Copyright (c) 2012, UT-Battelle, LLC\n"
"All rights reserved\n"
"[MPS++, Version 0.1]\n"
"\n"
"********************************************************************************\n"
"DISCLAIMER\n"
"\n"
"UT-Battelle, LLC AND THE GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL\n"
"WARRANTIES, BOTH EXPRESSED AND IMPLIED.  THERE ARE NO EXPRESS OR IMPLIED\n"
"WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE\n"
"USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR\n"
"OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED\n"
"RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.\n"
"THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS,\n"
"CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING\n"
"OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.\n"
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

typedef double RealType;
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
typedef PsimagLite::InputNg<Mpspp::InputCheck> InputNgType;
typedef PsimagLite::Geometry<RealType,Mpspp::ProgramGlobals> GeometryType;
typedef Mpspp::ParametersMpsSolver<RealType,InputNgType::Readable> ParametersMpsSolverType;

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

//	ParametersMpsSolverType mpsSolverParams(io);

//	ModelType model(mpsSolverParams,io,geometry,concurrency);

//	MatrixProductStateType psi; // initialize to something

//	MpsSolverType mpsSolver(mpsParams,model,concurrency);

//	mpsSolver.computeGroundState(psi);

//	const MatrixProductOperatorType& H = model.hamiltonian();

//	MatrixProductStateType hpsi = H*psi;

//	std::cout<<"Energy="<<scalarProduct(psi,hpsi)<<"\n";

//	std::cout<<"That's all folks!\n";
}

/*@}*/

