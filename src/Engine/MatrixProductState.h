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

#ifndef MATRIX_PRODUCT_STATE_H
#define MATRIX_PRODUCT_STATE_H

#include "ProgramGlobals.h"
#include "MpsFactor.h"

namespace Mpspp {

template<typename ComplexOrRealType_,typename SymmetryLocalType_>
class MatrixProductState {

	// FIXME: IDEA: PULL SYMMETRY OUT, PASS THROUGH FUNCTIONS

public:

	typedef SymmetryLocalType_ SymmetryLocalType;
	typedef typename SymmetryLocalType::IoInputType IoInputType;
	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef MpsFactor<ComplexOrRealType,SymmetryFactorType> MpsFactorType;
	typedef typename MpsFactorType::VectorType VectorType;
	typedef typename MpsFactorType::SparseMatrixType SparseMatrixType;

	MatrixProductState(size_t nsites)
	: nsites_(nsites),center_(0)
	{}

//	MatrixProductState(IoInputType& io)
//		: symmNonconst_(new SymmetryLocalType(io)),
//		  symm_(*symmNonconst_),
//		  nsites_(symm_(0).super().block().size()),
//		  center_(0)
//	{
//		io.rewind();
//		assert(nsites_>0);
//		for (size_t i=0;i<nsites_-1;i++) {
//			MpsFactorType f(io,symm_(i),i);
//			data_.push_back(f);
//		}
//	}

	~MatrixProductState()
	{
		for (size_t i=0;i<A_.size();i++)
			if (A_[i]) delete A_[i];
		for (size_t i=0;i<B_.size();i++)
			if (B_[i]) delete B_[i];
	}

	void growRight(size_t currentSite,const SymmetryLocalType& symm)
	{
		center_ = currentSite;
		MpsFactorType* mpsFactor = new MpsFactorType(currentSite,MpsFactorType::TYPE_B);
		mpsFactor->setRandom(currentSite,symm(currentSite));
		B_.push_back(mpsFactor);
		MpsFactorType* mpsFactor2 = new MpsFactorType(currentSite,MpsFactorType::TYPE_A);
		A_.push_back(mpsFactor2);
	}

	size_t center() const
	{
		return center_;
	}

	//! Returns the number of sites
	size_t sites() const
	{
		return nsites_;
	}

	//! tmpVec[i] --> M^\sigma2 _ {a1,a2}
	void update(size_t currentSite,const VectorType& v,size_t direction,size_t symmetrySector,const SymmetryFactorType& symm)
	{
		center_ = currentSite;
		if (direction==ProgramGlobals::TO_THE_RIGHT) {
			assert(currentSite<A_.size());
			A_[currentSite]->updateFromVector(v,symmetrySector,symm);
			std::cout<<"A_["<<currentSite<<"].row= "<<A_[currentSite]->operator()().row()<<"\n";
		} else {
			assert(currentSite<B_.size());
			B_[currentSite]->updateFromVector(v,symmetrySector,symm);
		}
	}

	const MpsFactorType& A(size_t site) const
	{
		assert(site<A_.size());
		return *(A_[site]);
	}

	const MpsFactorType& B(size_t site) const
	{
		assert(site<B_.size());
		return *(B_[site]);
	}


	typename ProgramGlobals::Real<ComplexOrRealType>::Type norm(size_t type) const
	{
		if (type==MpsFactorType::TYPE_B) {
			return normB();
		}
		return normA();
	}

	template<typename ComplexOrRealType2,typename SymmetryLocalType2>
	friend std::ostream& operator<<(std::ostream& os,const MatrixProductState<ComplexOrRealType2,SymmetryLocalType2>& mps);

private:

	typename ProgramGlobals::Real<ComplexOrRealType>::Type normA() const
	{
		throw std::runtime_error("normA(): unimplemented\n");
	}

	typename ProgramGlobals::Real<ComplexOrRealType>::Type normB() const
	{
		throw std::runtime_error("normB(): unimplemented\n");
	}

	// copy ctor:
	MatrixProductState(const MatrixProductState& other);

	// assignment
	MatrixProductState& operator=(const MatrixProductState& other);

//	const SymmetryLocalType& symm_;
	size_t nsites_;
	size_t center_;
	typename ProgramGlobals::Vector<MpsFactorType*>::Type B_;
	typename ProgramGlobals::Vector<MpsFactorType*>::Type A_;
}; // MatrixProductState

template<typename ComplexOrRealType,typename SymmetryLocalType>
std::ostream& operator<<(std::ostream& os,const MatrixProductState<ComplexOrRealType,SymmetryLocalType>& mps)
{
	os<<"nsites= "<<mps.nsites_<<" center="<<mps.center_;
	os<<"A_.size= "<<mps.A_.size()<<"\n";
	for (size_t i=0;i<mps.A_.size();i++)
		os<<*(mps.A_[i]);
	os<<"B_.size= "<<mps.B_.size()<<"\n";
	for (size_t i=0;i<mps.B_.size();i++)
		os<<*(mps.B_[i]);
	os<<"\n";
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // MATRIX_PRODUCT_STATE_H

