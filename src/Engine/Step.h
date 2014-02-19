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
#include "FermionSign.h"

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
	typedef typename MpsLocalType::VectorIntegerType VectorIntegerType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::SymmetryHelperType SymmetryHelperType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef typename ParametersSolverType::InputValidatorType InputValidatorType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	static int const MAX_ = 100;
	
public:

	Step(const ParametersSolverType& solverParams,
	     MpsLocalType& mps,
	     ContractedLocalType& contractedLocal,
	     const ModelType& model,
	     InputValidatorType& io)
	: progress_("Step"),
	  solverParams_(solverParams),
	  mps_(mps),
	  contractedLocal_(contractedLocal),
	  model_(model),
	  paramsForSolver_(io,"Lanczos"),
	  statePredictor_(),
	  truncation_(mps,contractedLocal,solverParams_.options.find("notruncation")==PsimagLite::String::npos)
	{}

	//! Moves the center of orthogonality by one to the left
	void moveLeft(SymmetryLocalType& symm,size_t currentSite,const FiniteLoop& finiteLoop)
	{
		if (currentSite==model_.geometry().numberOfSites()) return;
		VectorIntegerType quantumNumbers;
		model_.getOneSite(quantumNumbers,currentSite);

		symm.moveLeft(currentSite,quantumNumbers);
		if (currentSite==0) return;

		FermionSign<ModelType> fermionSign(model_,currentSite);
		SymmetryHelperType symmetryHelper(fermionSign,symm);
		internalmove(TO_THE_LEFT,symmetryHelper,currentSite);
		contractedLocal_.move(currentSite,TO_THE_LEFT,symmetryHelper);
		truncation_(symm,currentSite,ProgramGlobals::PART_RIGHT,finiteLoop.keptStates);
	}

	//! Moves the center of orthogonality by one to the right
	void moveRight(SymmetryLocalType& symm,size_t currentSite,const FiniteLoop& finiteLoop)
	{
		VectorIntegerType quantumNumbers;
		model_.getOneSite(quantumNumbers,currentSite);
		symm.moveRight(currentSite+1,quantumNumbers);
		if (currentSite+1==model_.geometry().numberOfSites()) return;

		FermionSign<ModelType> fermionSign(model_,currentSite);
		SymmetryHelperType symmetryHelper(fermionSign,symm);
		internalmove(TO_THE_RIGHT,symmetryHelper,currentSite+1);
		contractedLocal_.move(currentSite,TO_THE_RIGHT,symmetryHelper);
		truncation_(symm,currentSite,ProgramGlobals::PART_LEFT,finiteLoop.keptStates);
	}

	void grow(SymmetryLocalType& symm,size_t center)
	{
		size_t nsites = model_.geometry().numberOfSites();
		VectorIntegerType quantumNumbers;
		model_.getOneSite(quantumNumbers,center);
		symm.grow(center,quantumNumbers,nsites);
		mps_.grow(center,symm,quantumNumbers.size());

		FermionSign<ModelType> fermionSign(model_,center);
		SymmetryHelperType symmetryHelper(fermionSign,symm);
		contractedLocal_.grow(center,symmetryHelper,nsites);
		if (center==0) return;

		internalmove(TO_THE_RIGHT,symmetryHelper,center+1);

		FermionSign<ModelType> fermionSign2(model_,center+1);
		SymmetryHelperType symmetryHelper2(fermionSign2,symm);
		internalmove(TO_THE_LEFT,symmetryHelper2,center+1);

		truncation_(symm,center,ProgramGlobals::PART_LEFT,solverParams_.keptStatesInfinite);
		truncation_(symm,center+1,ProgramGlobals::PART_RIGHT,solverParams_.keptStatesInfinite);
	}

	void printReport(std::ostream& os) const
	{
		os<<"Nothing to report so far, except that I need a progress indicator\n";
	}

private:

	void internalmove(size_t direction,const SymmetryHelperType& symmetryHelper,size_t siteForSymm)
	{
		const SymmetryFactorType& symm = symmetryHelper.symmLocal()(siteForSymm);
		size_t currentSite = symmetryHelper.currentSite();
		size_t symmetrySector = getSymmetrySector(direction,symm.super());
		std::cerr<<"symmetrySector="<<symmetrySector<<"\n";
//		size_t total = symm.super().size();
		size_t total = symm.super().partitionSize(symmetrySector);
		VectorType v(total,0.0);
		RealType energy = internalmove(v,currentSite,direction,symmetryHelper,symmetrySector,siteForSymm);
		mps_.move(truncation_,currentSite,v,direction,symmetrySector,symm);
		statePredictor_.push(energy,v,symmetrySector);
		if (solverParams_.options.find("test")!=PsimagLite::String::npos)
			throw PsimagLite::LogicError
					 ("Exiting due to option test in the input file\n");
	}

	RealType internalmove(VectorType& tmpVec,size_t currentSite,size_t direction,const SymmetryHelperType& symmetryHelper,size_t symmetrySector,size_t siteForSymm)
	{
		const ParametersSolverType& solverParams = model_.solverParams();

		typedef InternalProductTemplate<typename VectorType::value_type,ModelType> InternalProductType;
		typedef PsimagLite::LanczosOrDavidsonBase<ParametersForSolverType,InternalProductType,VectorType> LanczosOrDavidsonBaseType;

		ReflectionSymmetryType *rs = 0;
		ModelHelperType modelHelper(contractedLocal_,symmetrySector,currentSite,direction,model_.hamiltonian()(currentSite),symmetryHelper,siteForSymm);
		InternalProductType lanczosHelper(&model_,&modelHelper,rs);

		LanczosOrDavidsonBaseType* lanczosOrDavidson = 0;

		bool useDavidson = (solverParams.options.find("useDavidson")!=PsimagLite::String::npos);
		if (useDavidson) {
			lanczosOrDavidson = new PsimagLite::DavidsonSolver<ParametersForSolverType,InternalProductType,VectorType>(lanczosHelper,paramsForSolver_);
		} else {
			lanczosOrDavidson = new PsimagLite::LanczosSolver<ParametersForSolverType,InternalProductType,VectorType>(lanczosHelper,paramsForSolver_);
		}

		RealType energyTmp = 0;
		assert(tmpVec.size()==lanczosHelper.rank());
		VectorType initialVector(lanczosHelper.rank());
		statePredictor_.createRandomVector(initialVector,0,initialVector.size());

		lanczosOrDavidson->computeGroundState(energyTmp,tmpVec,initialVector);
		if (lanczosOrDavidson) delete lanczosOrDavidson;
//		std::cout<<"Eigenstate\n";
//		std::cout<<tmpVec;
		return energyTmp;
	}

	size_t getSymmetrySector(size_t direction,const SymmetryComponentType& super) const
	{
		size_t sites = super.block().size();
		size_t targetQuantumNumber = getQuantumSector(sites,direction);
		size_t imin=0;
		size_t minDiff=0;
		for (size_t i=0;i<super.partitions()-1;i++) {
			size_t state = super.partitionOffset(i);
			size_t q = super.qn(state);
			size_t diff = (q<targetQuantumNumber) ? targetQuantumNumber - q : q-targetQuantumNumber;
			if (diff<minDiff || i==0) {
				imin = i;
				minDiff = diff;
			}
		}
		if (minDiff!=0) {
			std::cerr<<__FILE__<<" "<<__LINE__<<"\n";
			std::cerr<<"getSymmetrySector ";
			std::cerr<<"WARNING: minDiff="<<minDiff<<" is non zero\n";
			throw std::runtime_error("error\n");
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
		VectorIntegerType targetQuantumNumbers(solverParams_.targetQuantumNumbers.size());
		for (size_t ii=0;ii<targetQuantumNumbers.size();ii++) 
			targetQuantumNumbers[ii]=size_t(round(solverParams_.targetQuantumNumbers[ii]*sites));
		return getQuantumSector(targetQuantumNumbers,direction);
	}

	size_t getQuantumSectorUd(size_t sites,size_t direction) const
	{
		VectorIntegerType targetQuantumNumbers(2);

		size_t nsites = model_.geometry().numberOfSites();
		targetQuantumNumbers[0]=size_t(static_cast<RealType>(solverParams_.electronsUp*sites)/nsites);
		targetQuantumNumbers[1]=size_t(static_cast<RealType>(solverParams_.electronsDown*sites)/nsites);

		return getQuantumSector(targetQuantumNumbers,direction);
	}

	size_t getQuantumSector(const VectorIntegerType& targetQuantumNumbers,size_t direction) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Integer target quantum numbers are: ";
		for (size_t ii=0;ii<targetQuantumNumbers.size();ii++)
			msg<<targetQuantumNumbers[ii]<<" ";
		progress_.printline(msg,std::cout);
		//if (direction==INFINITE) io_.printVector(targetQuantumNumbers,"TargetedQuantumNumbers");
		return encodeQuantumNumber(targetQuantumNumbers);
	}
	
	static size_t encodeQuantumNumber(const VectorIntegerType& v)
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
	ParametersForSolverType paramsForSolver_;
	StatePredictorType statePredictor_;
	TruncationType truncation_;
}; // Step

} // namespace Mpspp

/*@}*/
#endif // STEP_H

