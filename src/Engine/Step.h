/*
Copyright (c) 2012, UT-Battelle, LLC
All rights reserved

[MPS++, Version 0.1]
[by K. Al-Hassanieh, Oak Ridge National Laboratory]
[by J. Rincon, Oak Ridge National Laboratory]
[by G.A., Oak Ridge National Laboratory]

See full open source LICENSE under file LICENSE
in the root directory of this distribution.

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup MPSPP */
/*@{*/

#ifndef STEP_H
#define STEP_H

#include "ContractedPart.h"
#include "ProgressIndicator.h"
#include "ParametersForSolver.h"
#include "LanczosOrDavidsonBase.h"
#include "LanczosSolver.h"
#include "DavidsonSolver.h"
#include "LeftRightSuper.h"
#include <vector>

namespace Mpspp {

template<typename ModelType,
		 template<typename,typename> class InternalProductTemplate>
class Step {

	typedef typename ModelType::MatrixProductOperatorType MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef typename ModelType::ParametersSolverType ParametersSolverType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ModelType::VectorType VectorType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

public:

	typedef LeftRightSuper<MatrixProductOperatorType,VectorType,RealType> LeftRightSuperType;

	Step(LeftRightSuperType& lrs,const ModelType& model)
	: progress_("Step",0),
	  lrs_(lrs),
	  model_(model)
	{}

	//! Moves the center of orthogonality by one to the right
	void moveRight()
	{
		internalUpdate(TO_THE_RIGHT); // <--  From cL and cR construct a new A, only A changes here
		lrs_.contractedLeft().update(lrs_.A(),TO_THE_RIGHT);
	}

	//! Moves the center of orthogonality by one to the left
	void moveLeft()
	{
		internalUpdate(TO_THE_LEFT); // <-- From cL and cR construct a new B, only B changes here
		lrs_.contractedRight().update(lrs_.B(),TO_THE_LEFT);
	}

	void printReport(std::ostream& os) const
	{
		os<<"Nothing to report so far, except that I need a progress indicator\n";
	}

private:

	void internalUpdate(size_t direction)
	{
		const ParametersSolverType& solverParams = model_.solverParams();

		typedef InternalProductTemplate<typename VectorType::value_type,ModelType> InternalProductType;
		typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
		typedef PsimagLite::LanczosOrDavidsonBase<ParametersForSolverType,InternalProductType,VectorType> LanczosOrDavidsonBaseType;
		typedef typename ModelType::ModelHelperType ModelHelperType;

		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to set symmetry sector here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());

		ReflectionSymmetryType *rs = 0;
		ModelHelperType modelHelper(direction);
		typename LanczosOrDavidsonBaseType::MatrixType lanczosHelper(&model_,&modelHelper,rs);

		RealType eps=ProgramGlobals::LanczosTolerance;
		int iter=ProgramGlobals::LanczosSteps;

		ParametersForSolverType params;
		params.steps = iter;
		params.tolerance = eps;
		params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;
		params.options= solverParams.options;
		params.lotaMemory=false; //!(parameters_.options.find("DoNotSaveLanczosVectors")!=std::string::npos);

		LanczosOrDavidsonBaseType* lanczosOrDavidson = 0;

		bool useDavidson = (solverParams.options.find("useDavidson")!=std::string::npos);
		if (useDavidson) {
			lanczosOrDavidson = new PsimagLite::DavidsonSolver<ParametersForSolverType,InternalProductType,VectorType>(lanczosHelper,params);
		} else {
			lanczosOrDavidson = new PsimagLite::LanczosSolver<ParametersForSolverType,InternalProductType,VectorType>(lanczosHelper,params);
		}

		RealType energyTmp = 0;
		VectorType tmpVec(lanczosHelper.rank());
		VectorType initialVector(lanczosHelper.rank());
		std::string str2(__FILE__);
		str2 += " " + ttos(__LINE__) + "\n";
		str2 += "Initial vector must be set here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str2.c_str());

		lanczosOrDavidson->computeGroundState(energyTmp,tmpVec,initialVector);
		if (lanczosOrDavidson) delete lanczosOrDavidson;

		vector2Mps((direction==TO_THE_RIGHT) ? lrs_.A() : lrs_.B(),tmpVec);
	}

	//! tmpVec[i] --> M^\sigma2 _ {a1,a2}
	void vector2Mps(MatrixProductStateType& AorB,const VectorType& tmpVec)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need vector2Mps(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	PsimagLite::ProgressIndicator progress_;
	LeftRightSuperType& lrs_;
	const ModelType& model_;
}; // Step

} // namespace Mpspp

/*@}*/
#endif // STEP_H

