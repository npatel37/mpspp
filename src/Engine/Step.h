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
#include <vector>
#include "StatePredictor.h"

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
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::SymmetryComponentType SymmetryComponentType;
	typedef StatePredictor<RealType> StatePredictorType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	static int const MAX_ = 100;
	
public:

	Step(const ParametersSolverType& solverParams,LeftRightSuperType& lrs,const ModelType& model)
	: progress_("Step",0),
	  solverParams_(solverParams),
	  lrs_(lrs),
	  model_(model),
	  statePredictor_()
	{}

	//! Moves the center of orthogonality by one to the right
	void moveRight(size_t currentSite)
	{
		internalUpdate(currentSite,TO_THE_RIGHT); // <--  From cL and cR construct a new A, only A changes here
		lrs_.updateContracted(currentSite,lrs_.abState(),TO_THE_RIGHT);
	}

	//! Moves the center of orthogonality by one to the left
	void moveLeft(SymmetryLocalType& symm,size_t currentSite)
	{
		std::vector<size_t> quantumNumbers;
		size_t hilbert = 0;
		model_.getOneSite(hilbert,quantumNumbers,currentSite);
		symm.moveLeft(hilbert,currentSite,quantumNumbers);
		lrs_.moveLeft(currentSite); // computes A, computes L
		internalUpdate(currentSite,TO_THE_LEFT); // <-- From cL and cR construct a new B, only B changes here
		lrs_.updateContracted(currentSite,TO_THE_LEFT);
	}

	void growRight(SymmetryLocalType& symm,size_t currentSite)
	{
		std::vector<size_t> quantumNumbers;
		size_t hilbert = 0;
		model_.getOneSite(hilbert,quantumNumbers,currentSite);
		symm.growRight(hilbert,currentSite,quantumNumbers); // grows symm
		lrs_.growRight(currentSite); // grows B, computes R
		internalUpdate(currentSite,TO_THE_RIGHT); // <--  From cL and cR construct a new A, only A changes here
		lrs_.updateContracted(currentSite,TO_THE_RIGHT);
	}

	void printReport(std::ostream& os) const
	{
		os<<"Nothing to report so far, except that I need a progress indicator\n";
	}

private:

	void internalUpdate(size_t currentSite,size_t direction)
	{
		const ParametersSolverType& solverParams = model_.solverParams();

		typedef InternalProductTemplate<typename VectorType::value_type,ModelType> InternalProductType;
		typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
		typedef PsimagLite::LanczosOrDavidsonBase<ParametersForSolverType,InternalProductType,VectorType> LanczosOrDavidsonBaseType;
		typedef typename ModelType::ModelHelperType ModelHelperType;

		size_t symmetrySector = getSymmetrySector(direction);

		ReflectionSymmetryType *rs = 0;
		ModelHelperType modelHelper(lrs_,symmetrySector,currentSite,direction,model_.hamiltonian()(currentSite));
		InternalProductType lanczosHelper(&model_,&modelHelper,rs);

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
		statePredictor_.createRandomVector(initialVector,0,initialVector.size());

		lanczosOrDavidson->computeGroundState(energyTmp,tmpVec,initialVector);
		if (lanczosOrDavidson) delete lanczosOrDavidson;

		lrs_.updateMps(currentSite,tmpVec,direction,symmetrySector);
	}

	size_t getSymmetrySector(size_t direction) const
	{
		size_t center = lrs_.abState().center();
		size_t sites = lrs_.symmetry()(center).super().block().size();
		size_t targetQuantumNumber = getQuantumSector(sites,direction);
		
		SymmetryComponentType super = lrs_.symmetry()(center).super();
		for (size_t i=0;i<super.partitions()-1;i++) {
			size_t state = super.partitionOffset(i);
			size_t q = super.qn(state);
			if (q == targetQuantumNumber) return i;
		}
		assert(false);
		throw std::runtime_error("getSymmetrySector\n");
		return -1;
	}
	
	size_t getQuantumSector(size_t sites,size_t direction) const
	{
		return (solverParams_.targetQuantumNumbers.size()>0) ?
				getQuantumSectorT(sites,direction) :
				getQuantumSectorUd(sites,direction);
	}

	size_t getQuantumSectorT(size_t sites,size_t direction) const
	{
		std::vector<size_t> targetQuantumNumbers(solverParams_.targetQuantumNumbers.size());
		for (size_t ii=0;ii<targetQuantumNumbers.size();ii++) 
			targetQuantumNumbers[ii]=size_t(round(solverParams_.targetQuantumNumbers[ii]*sites));
		return getQuantumSector(targetQuantumNumbers,direction);
	}

	size_t getQuantumSectorUd(size_t sites,size_t direction) const
	{
		std::vector<size_t> targetQuantumNumbers(2);

		size_t nsites = model_.geometry().numberOfSites();
		targetQuantumNumbers[0]=static_cast<RealType>(solverParams_.electronsUp*sites)/nsites;
		targetQuantumNumbers[1]=static_cast<RealType>(solverParams_.electronsDown*sites)/nsites;

		return getQuantumSector(targetQuantumNumbers,direction);
	}

	size_t getQuantumSector(const std::vector<size_t>& targetQuantumNumbers,size_t direction) const
	{
		std::ostringstream msg;
		msg<<"Integer target quantum numbers are: ";
		for (size_t ii=0;ii<targetQuantumNumbers.size();ii++)
			msg<<targetQuantumNumbers[ii]<<" ";
		progress_.printline(msg,std::cout);
		//if (direction==INFINITE) io_.printVector(targetQuantumNumbers,"TargetedQuantumNumbers");
		return encodeQuantumNumber(targetQuantumNumbers);
	}
	
	static size_t encodeQuantumNumber(const std::vector<size_t>& v)
	{
		size_t x= v[0] + v[1]*MAX_;
		if (v.size()==3) x += v[2]*MAX_*MAX_;
		return x;
	}

	PsimagLite::ProgressIndicator progress_;
	const ParametersSolverType& solverParams_;
	LeftRightSuperType& lrs_;
	const ModelType& model_;
	StatePredictorType statePredictor_;
}; // Step

} // namespace Mpspp

/*@}*/
#endif // STEP_H

