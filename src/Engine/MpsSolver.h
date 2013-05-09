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
	typedef typename ModelBaseType::MpoLocalType MpoLocalType;
	typedef typename MpoLocalType::MpsLocalType MpsLocalType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ParametersSolverType::FiniteLoopsType FiniteLoopsType;
	typedef typename MpsLocalType::VectorIntegerType VectorIntegerType;
	typedef typename ModelBaseType::ContractedLocalType ContractedLocalType;
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
		MpsLocalType psi(model_.geometry().numberOfSites());
		ContractedLocalType contracted(psi,model_.hamiltonian());

		size_t center = 0;

		growLattice(psi,contracted,center,symm);

		if (solverParams_.options.find("nofiniteloops")!=PsimagLite::String::npos)
			return;

		StepType step(solverParams_,psi,contracted,model_);
		finiteLoops(step,center,symm);
	}

private:

	void growLattice(MpsLocalType& psi,ContractedLocalType& contracted,size_t& center,SymmetryLocalType& symm)
	{
		size_t nsites = model_.geometry().numberOfSites();

		StepType step(solverParams_,psi,contracted,model_);

		if (nsites&1) {
			throw std::runtime_error("growLattice: nsites must be even\n");
		}
		size_t total = size_t(nsites/2);
		for (size_t i=0;i<total;i++) {
			center=i;
			step.grow(symm,center);
		}
	}

	void finiteLoops(StepType& step,size_t& center,SymmetryLocalType& symm)
	{
		size_t nsites = model_.geometry().numberOfSites();

		for (size_t finiteLoop=0;finiteLoop<solverParams_.finiteLoops.size();finiteLoop++) {
			FiniteLoop fl = solverParams_.finiteLoops[finiteLoop];
			size_t counter = abs(fl.stepLength);
			if (finiteLoop==0) {
				if (fl.stepLength>0) {
					center++;
				}
			}
			while(true) {
				if (fl.stepLength<0) {
					step.moveLeft(symm,center,fl);
					printProgress(symm,center,"to the left");
					if (center==0) break;
					if (counter==0) break;
					counter--;
					if (counter==0) break;
					center--;
				} else {
					step.moveRight(symm,center,fl);
					printProgress(symm,center,"to the right");
					center++;
					if (center==nsites) break;
					if (counter==0) break;
					counter--;
					if (counter==0) break;
				}
			}
		}
	}

	void printProgress(const SymmetryLocalType& symm,size_t center,const PsimagLite::String& direction) const
	{
		PsimagLite::OstringStream msg;
		msg<<"center="<<center<<" direction="<<direction;
		msg<<" left space= "<<symm(center).left().size();
		msg<<" right space="<<symm(center).right().size()<<"\n";
		progress_.printline(msg,std::cout);
	}

	void printMemoryUsage() const
	{
		PsimagLite::MemoryUsage musage;
		PsimagLite::String vmPeak = musage.findEntry("VmPeak:");
		PsimagLite::String vmSize = musage.findEntry("VmSize:");
		PsimagLite::OstringStream msg;
		msg<<"Current virtual memory is "<<vmSize<<" maximum was "<<vmPeak;
		progress_.printline(msg,std::cout);
		PsimagLite::OstringStream msg2;
		msg2<<"Amount of time scheduled (user plus system): "<<musage.time()<<" clock ticks";
		progress_.printline(msg2,std::cout);
	}

	const ParametersSolverType& solverParams_;
	const ModelBaseType& model_;
	ConcurrencyType& concurrency_;
	PsimagLite::ProgressIndicator progress_;
	int stepCurrent_;
	VectorIntegerType sitesIndices_;
}; // MpsSolver

} // namespace Mpspp

/*@}*/
#endif // MPS_SOLVER_H

