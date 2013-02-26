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

#ifndef CONTRACTED_LOCAL_H
#define CONTRACTED_LOCAL_H

#include "ProgramGlobals.h"
#include "ContractedFactor.h"

namespace Mpspp {

template<typename MatrixProductOperatorType_>
class ContractedLocal {

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT, PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:

	typedef MatrixProductOperatorType_ MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MatrixProductOperatorType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::SymmetryLocalType SymmetryLocalType;
	typedef typename MpsLocalType::ComplexOrRealType ComplexOrRealType;
	typedef typename ProgramGlobals::CrsMatrix<ComplexOrRealType>::Type SparseMatrixType;
	typedef ContractedFactor<MatrixProductOperatorType> ContractedFactorType;

	ContractedLocal(const MpsLocalType& abState,const MatrixProductOperatorType& h)
		: abState_(abState),
		  h_(h),
		  R_(abState.sites(),ProgramGlobals::PART_RIGHT),
		  L_(abState.sites(),ProgramGlobals::PART_LEFT)
	{}

	void grow(size_t currentSite,const SymmetryLocalType& symm,size_t nsites)
	{
		L_[currentSite+1].build(abState_.A(currentSite),h_(currentSite),L_[currentSite],symm(currentSite+1));

		R_[currentSite+1].build(abState_.B(currentSite),h_(nsites-1-currentSite),R_[currentSite],symm(currentSite+1));
	}

	//! From As (or Bs) and Ws reconstruct *this
	void move(size_t currentSite,size_t direction,const SymmetryLocalType& symm)
	{
		if (direction==TO_THE_RIGHT) {
			moveLeft(currentSite,abState_,symm);
		} else {
			moveRight(currentSite,abState_,symm);
		}
	}

	void truncate(size_t site,size_t part,size_t cutoff,size_t nsites)
	{
		if (part==ProgramGlobals::PART_LEFT) {
			if (site+1>=L_.size()) return;
			L_[site+1].truncate(cutoff);
		} else {
			size_t siteToSet = nsites - site;
			R_[siteToSet].truncate(cutoff);
		}
	}

	const ContractedFactorType& operator()(size_t currentSite,size_t leftOrRight) const
	{
		if (leftOrRight == PART_LEFT) {
			assert(currentSite<L_.size());
		} else {
			assert(currentSite<R_.size());
		}

		return (leftOrRight == PART_LEFT) ? L_[currentSite] : R_[currentSite];
	}

	template<typename MatrixProductOperatorType2>
	friend std::ostream& operator<<(std::ostream& os,const ContractedLocal<MatrixProductOperatorType2>& contractedLocal);

private:

	void moveLeft(size_t currentSite,const MpsLocalType& abState,const SymmetryLocalType& symm)
	{
		assert(currentSite+1<L_.size());
		L_[currentSite+1].move(abState.A(currentSite),h_(currentSite),L_[currentSite],symm(currentSite+1));
		std::cout<<"set L_["<<(currentSite+1)<<"]="<<L_[currentSite+1].row()<<"\n";
	}

	void moveRight(size_t currentSite,const MpsLocalType& abState,const SymmetryLocalType& symm)
	{
		size_t nsites = symm(currentSite).super().block().size();
		size_t siteToSet = nsites - currentSite;
		if (siteToSet==R_.size()) return;
		assert(siteToSet<R_.size());
		R_[siteToSet].move(abState.B(siteToSet-1),h_(currentSite),R_[siteToSet-1],symm(currentSite));
		std::cout<<"set R_["<<siteToSet<<"]="<<R_[siteToSet].row()<<"\n";
	}

	const MpsLocalType& abState_;
	const MatrixProductOperatorType& h_;
	typename ProgramGlobals::Vector<ContractedFactorType>::Type R_;
	typename ProgramGlobals::Vector<ContractedFactorType>::Type L_;

}; // ContractedLocal

template<typename MatrixProductOperatorType>
std::ostream& operator<<(std::ostream& os,const ContractedLocal<MatrixProductOperatorType>& contractedLocal)
{
	os<<"ContractedLocal: right size="<<contractedLocal.R_.size()<<"\n";
	for (size_t i=0;i<contractedLocal.R_.size();i++)
		os<<contractedLocal.R_[i];
	os<<"ContractedLocal: left size="<<contractedLocal.L_.size()<<"\n";
	for (size_t i=0;i<contractedLocal.L_.size();i++)
		os<<contractedLocal.L_[i];
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_LOCAL_H

