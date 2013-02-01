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

#ifndef MPS_SOLVER_H
#define MPS_SOLVER_H

#include "Step.h"
#include "ProgramGlobals.h"
#include "ProgressIndicator.h"
#include "MemoryUsage.h"

namespace Mpspp {

template<typename ModelBaseType,
template<typename,typename> class InternalProductTemplate>
class MpsSolver {

	typedef typename ModelBaseType::ParametersSolverType ParametersSolverType;
	typedef typename ModelBaseType::InputValidatorType InputValidatorType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelBaseType::ConcurrencyType ConcurrencyType;
	typedef typename ModelBaseType::MatrixProductOperatorType MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ParametersSolverType::FiniteLoopsType FiniteLoopsType;
//	typedef typename ModelBaseType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelBaseType::ContractedPartType ContractedPartType;
	typedef typename ModelBaseType::SymmetryLocalType SymmetryLocalType;
	typedef Step<ModelBaseType,InternalProductTemplate> StepType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

public:

	MpsSolver(const ParametersSolverType& solverParams,
			  const ModelBaseType& model,
			  ConcurrencyType& concurrency)
		: solverParams_(solverParams),
		  model_(model),
		  concurrency_(concurrency),
		  progress_("MpsSolver",concurrency.rank()),
		  stepCurrent_(0),
		  sitesIndices_(model.geometry().numberOfSites())
	{
		for (size_t i=0;i<sitesIndices_.size();i++) sitesIndices_[i] = i;
	}

	void computeGroundState()
	{
		SymmetryLocalType symm;
		MatrixProductStateType psi(model_.geometry().numberOfSites());
		ContractedPartType contracted(psi,model_.hamiltonian());

		size_t center = 0;

		growLattice(psi,contracted,center,symm);

		StepType step(solverParams_,psi,contracted,model_);
		finiteLoops(step,center,symm);
	}

private:

	void growLattice(MatrixProductStateType& psi,ContractedPartType& contracted,size_t& center,SymmetryLocalType& symm)
	{
		size_t nsites = model_.geometry().numberOfSites();
		for (size_t i=0;i<nsites;i++) {
			center=i;
			std::vector<size_t> quantumNumbers;
			model_.getOneSite(quantumNumbers,center);
			symm.growRight(center,quantumNumbers,nsites); // grows symm
		}

		std::cout<<symm;

		for (size_t i=0;i<nsites;i++) {
			center=i;
			psi.growRight(center,symm);
		}

		std::cout<<psi;

//		RealType tmp = psi.norm(MatrixProductStateType::MpsFactorType::TYPE_B);
//		std::cout<<tmp<<"\n";

		for (size_t i=0;i<nsites;i++) {
			center = nsites-i-1;
			contracted.growRight(center,symm,nsites);
		}
		std::cout<<contracted;
	}

	void finiteLoops(StepType& step,size_t& center,SymmetryLocalType& symm)
	{
		size_t nsites = model_.geometry().numberOfSites();

		for (size_t finiteLoop=0;finiteLoop<solverParams_.finiteLoops.size();finiteLoop++) {
			FiniteLoop fl = solverParams_.finiteLoops[finiteLoop];
			size_t counter = abs(fl.stepLength);
			while(true) {
				if (fl.stepLength<0) {
					step.moveLeft(symm,center);
					printProgress(symm,center);
					if (center==0) break;
					if (counter==0) break;
					center--;
				} else {
					step.moveRight(symm,center);
					printProgress(symm,center);
					if (center+2==nsites) break;
					center++;
					if (counter==0) break;
					counter--;
				}
			}
		}
	}

	void printProgress(const SymmetryLocalType& symm,size_t center) const
	{
		std::ostringstream msg;
		msg<<"center="<<center<<" ";
		msg<<" left space= "<<symm(center).left().size();
		msg<<" right space="<<symm(center).right().size()<<"\n";
		progress_.printline(msg,std::cout);
	}

	void printMemoryUsage() const
	{
		PsimagLite::MemoryUsage musage;
		std::string vmPeak = musage.findEntry("VmPeak:");
		std::string vmSize = musage.findEntry("VmSize:");
		std::ostringstream msg;
		msg<<"Current virtual memory is "<<vmSize<<" maximum was "<<vmPeak;
		progress_.printline(msg,std::cout);
		std::ostringstream msg2;
		msg2<<"Amount of time scheduled (user plus system): "<<musage.time()<<" clock ticks";
		progress_.printline(msg2,std::cout);
	}

	const ParametersSolverType& solverParams_;
	const ModelBaseType& model_;
	ConcurrencyType& concurrency_;
	PsimagLite::ProgressIndicator progress_;
	int stepCurrent_;
	std::vector<size_t> sitesIndices_;
}; // MpsSolver

} // namespace Mpspp

/*@}*/
#endif // MPS_SOLVER_H

