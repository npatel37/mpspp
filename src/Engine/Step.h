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

//#include "ContractedLocal.h"
#include "ProgressIndicator.h"
#include "ParametersForSolver.h"
#include "LanczosOrDavidsonBase.h"
#include "LanczosSolver.h"
#include "DavidsonSolver.h"
#include <vector>
#include "StatePredictor.h"
#include "Truncation.h"

namespace Mpspp {

template<typename ModelType,
         template<typename,typename> class InternalProductTemplate>
class Step {

	typedef typename ModelType::MpoLocalType MpoLocalType;
	typedef typename MpoLocalType::MpsLocalType MpsLocalType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef typename ModelType::ParametersSolverType ParametersSolverType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ModelType::VectorType VectorType;
	typedef typename ModelType::ContractedLocalType ContractedLocalType;
	typedef typename ContractedLocalType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::SymmetryComponentType SymmetryComponentType;
	typedef StatePredictor<RealType,VectorType> StatePredictorType;
	typedef typename MpsLocalType::VectorRealType VectorRealType;
	typedef Truncation<ContractedLocalType> TruncationType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	static int const MAX_ = 100;
	
public:

	Step(const ParametersSolverType& solverParams,
		 MpsLocalType& mps,
		 ContractedLocalType& contractedLocal,
		 const ModelType& model)
	: progress_("Step",0),
	  solverParams_(solverParams),
	  mps_(mps),
	  contractedLocal_(contractedLocal),
	  model_(model),
	  statePredictor_(),
	  truncation_(mps,contractedLocal)
	{}

	//! Moves the center of orthogonality by one to the left
	void moveLeft(SymmetryLocalType& symm,size_t currentSite,const FiniteLoop& finiteLoop)
	{
		if (currentSite+1==model_.geometry().numberOfSites()) return;
		std::vector<size_t> quantumNumbers;
		model_.getOneSite(quantumNumbers,currentSite);

		symm.moveLeft(currentSite,quantumNumbers);
//		TruncationType truncation(symm(currentSite).left().size());
		internalmove(currentSite,TO_THE_LEFT,symm(currentSite));
		contractedLocal_.move(currentSite,TO_THE_LEFT,symm);
//		truncation.truncate(symm(currentSite).right().size());
		truncation_(symm,currentSite,ProgramGlobals::PART_RIGHT,finiteLoop.keptStates);
	}

	//! Moves the center of orthogonality by one to the right
	void moveRight(SymmetryLocalType& symm,size_t currentSite,const FiniteLoop& finiteLoop)
	{
		std::vector<size_t> quantumNumbers;
		model_.getOneSite(quantumNumbers,currentSite);

		symm.moveRight(currentSite+1,quantumNumbers);
//		std::cout<<"normB="<<mps_.norm(MpsLocalType::MpsFactorType::TYPE_B,symm)<<" ";
//		std::cout<<"normA="<<mps_.norm(MpsLocalType::MpsFactorType::TYPE_A,symm)<<"\n";
//		TruncationType truncation(symm(currentSite).right().size());
		internalmove(currentSite,TO_THE_RIGHT,symm(currentSite+1));
		contractedLocal_.move(currentSite,TO_THE_RIGHT,symm);
//		truncation.truncate(symm(currentSite).left().size());
		truncation_(symm,currentSite,ProgramGlobals::PART_LEFT,finiteLoop.keptStates);
	}

	void grow(SymmetryLocalType& symm,size_t center)
	{
		size_t nsites = model_.geometry().numberOfSites();
		std::vector<size_t> quantumNumbers;
		model_.getOneSite(quantumNumbers,center);
		symm.grow(center,quantumNumbers,nsites);
		mps_.grow(center,symm,quantumNumbers.size());
		contractedLocal_.grow(center,symm,nsites);
		if (center==0) return;
//		TruncationType truncation(symm(center).right().size());
		internalmove(center,TO_THE_RIGHT,symm(center+1));
		internalmove(center+1,TO_THE_LEFT,symm(center+1));
//		truncation_(symm,center,ProgramGlobals::PART_LEFT,solverParams_.keptStatesInfinite);
	}

	void printReport(std::ostream& os) const
	{
		os<<"Nothing to report so far, except that I need a progress indicator\n";
	}

private:

	void internalmove(size_t currentSite,size_t direction,const SymmetryFactorType& symm)
	{
		size_t symmetrySector = getSymmetrySector(direction,symm.super());
		std::cerr<<"symmetrySector="<<symmetrySector<<"\n";
//		size_t total = symm.super().size();
		size_t total = symm.super().partitionSize(symmetrySector);
		VectorType v(total,0.0);
		RealType energy = internalmove(v,currentSite,direction,symm,symmetrySector);
		mps_.move(truncation_,currentSite,v,direction,symmetrySector,symm);
		statePredictor_.push(energy,v,symmetrySector);
	}

	RealType internalmove(VectorType& tmpVec,size_t currentSite,size_t direction,const SymmetryFactorType& symm,size_t symmetrySector)
	{
		const ParametersSolverType& solverParams = model_.solverParams();

		typedef InternalProductTemplate<typename VectorType::value_type,ModelType> InternalProductType;
		typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
		typedef PsimagLite::LanczosOrDavidsonBase<ParametersForSolverType,InternalProductType,VectorType> LanczosOrDavidsonBaseType;
		typedef typename ModelType::ModelHelperType ModelHelperType;

		ReflectionSymmetryType *rs = 0;
		size_t hamiltonianSite = (direction == TO_THE_RIGHT) ? currentSite : currentSite;
		ModelHelperType modelHelper(contractedLocal_,symmetrySector,currentSite,direction,model_.hamiltonian()(hamiltonianSite),symm);
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
		assert(tmpVec.size()==lanczosHelper.rank());
		VectorType initialVector(lanczosHelper.rank());
		statePredictor_.createRandomVector(initialVector,0,initialVector.size());

		lanczosOrDavidson->computeGroundState(energyTmp,tmpVec,initialVector);
		if (lanczosOrDavidson) delete lanczosOrDavidson;
		std::cout<<"Eigenstate\n";
		std::cout<<tmpVec;
		return energyTmp;
	}

	size_t getSymmetrySector(size_t direction,const SymmetryComponentType& super) const
	{
		size_t sites = super.block().size();
		size_t targetQuantumNumber = getQuantumSector(sites,direction);
		size_t imin=0;
		size_t minDiff=1e10;
		for (size_t i=0;i<super.partitions()-1;i++) {
			size_t state = super.partitionOffset(i);
			size_t q = super.qn(state);
			size_t diff = (q<targetQuantumNumber) ? targetQuantumNumber - q : q-targetQuantumNumber;
			if (diff<minDiff) {
				imin = i;
				minDiff = diff;
			}
		}
		if (minDiff!=0) {
			std::cerr<<__FILE__<<" "<<__LINE__<<"\n";
			std::cerr<<"getSymmetrySector ";
			std::cerr<<"WARNING: minDiff="<<minDiff<<" is non zero\n";
		}
		return imin;
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
	MpsLocalType& mps_;
	ContractedLocalType& contractedLocal_;
	const ModelType& model_;
	StatePredictorType statePredictor_;
	TruncationType truncation_;
}; // Step

} // namespace Mpspp

/*@}*/
#endif // STEP_H

