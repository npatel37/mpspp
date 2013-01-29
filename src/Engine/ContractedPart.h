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

#ifndef CONTRACTED_PART_H
#define CONTRACTED_PART_H

#include "ProgramGlobals.h"
#include "ContractedFactor.h"

namespace Mpspp {

template<typename MatrixProductOperatorType>
class ContractedPart {

	typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;
	typedef typename MatrixProductStateType::SymmetryLocalType SymmetryLocalType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MatrixProductStateType::ComplexOrRealType ComplexOrRealType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT, PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:

	typedef typename ProgramGlobals::CrsMatrix<ComplexOrRealType>::Type SparseMatrixType;
	typedef ContractedFactor<MatrixProductOperatorType> ContractedFactorType;

	ContractedPart(const MatrixProductStateType& abState,const MatrixProductOperatorType& h)
		: abState_(abState),h_(h)
	{}

	void growRight(size_t currentSite,const SymmetryLocalType& symm,size_t nsites)
	{
		ContractedFactorType cf(abState_.B(currentSite),h_(currentSite),currentSite,ProgramGlobals::PART_RIGHT,0,symm(currentSite));
		R_.push_back(cf);

		if (currentSite+2==nsites) {
			ContractedFactorType cf3(currentSite,ProgramGlobals::PART_RIGHT);
			R_.push_back(cf3);
		}

		if (currentSite>0) return;

		ContractedFactorType cf2(currentSite,ProgramGlobals::PART_LEFT);
		L_.push_back(cf2);
	}

	//! From As (or Bs) and Ws reconstruct *this
	void update(size_t currentSite,const MatrixProductStateType& abState,size_t direction,const SymmetryLocalType& symm)
	{
		if (direction==TO_THE_RIGHT) {
			updateLeft(currentSite,abState,symm);
		} else {
			updateRight(currentSite,abState,symm);
		}
	}

	const ContractedFactorType& operator()(size_t currentSite,size_t leftOrRight) const
	{
		if (leftOrRight == PART_LEFT) {
			assert(currentSite<L_.size());
//			std::cout<<L_;
		} else {
			assert(currentSite<R_.size());
		}

		return (leftOrRight == PART_LEFT) ? L_[currentSite] : R_[currentSite];
	}

	template<typename MatrixProductOperatorType_>
	friend std::ostream& operator<<(std::ostream& os,const ContractedPart<MatrixProductOperatorType_>& contractedPart);

private:

	void updateLeft(size_t currentSite,const MatrixProductStateType& abState,const SymmetryLocalType& symm)
	{
		assert(L_.size()>0);
		if (currentSite+1>=L_.size()) {
			ContractedFactorType* dataPrev = &(L_[L_.size()-1]);
			ContractedFactorType cf(abState.A(currentSite),h_(currentSite),currentSite,ProgramGlobals::PART_LEFT,dataPrev,symm(currentSite));
			L_.push_back(cf);
			return;
		}
		assert(currentSite<L_.size());
		assert(currentSite>0);
		L_[currentSite].update(abState.A(currentSite),h_(currentSite),L_[currentSite-1],symm(currentSite));
	}

	void updateRight(size_t currentSite,const MatrixProductStateType& abState,const SymmetryLocalType& symm)
	{
		assert(currentSite+1<R_.size());
		R_[currentSite].update(abState.B(currentSite),h_(currentSite),R_[currentSite+1],symm(currentSite));
	}

	const MatrixProductStateType& abState_;
	const MatrixProductOperatorType& h_;
	typename ProgramGlobals::Vector<ContractedFactorType>::Type R_;
	typename ProgramGlobals::Vector<ContractedFactorType>::Type L_;

}; // ContractedPart

template<typename MatrixProductOperatorType>
std::ostream& operator<<(std::ostream& os,const ContractedPart<MatrixProductOperatorType>& contractedPart)
{
	os<<"ContractedPart: right size="<<contractedPart.R_.size()<<"\n";
	for (size_t i=0;i<contractedPart.R_.size();i++)
		os<<contractedPart.R_[i];
	os<<"ContractedPart: left size="<<contractedPart.L_.size()<<"\n";
	for (size_t i=0;i<contractedPart.L_.size();i++)
		os<<contractedPart.L_[i];
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_PART_H

