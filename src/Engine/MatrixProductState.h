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

	MatrixProductState(size_t nsites,const SymmetryLocalType& symm)
	: nsites_(nsites),symmNonconst_(0),symm_(symm),center_(0)
	{}

//	MatrixProductState(size_t nsites)
//	: nsites_(nsites),symmNonconst_(new SymmetryLocalType()),symm_(*symmNonconst_),center_(0)
//	{
//		std::string str(__FILE__);
//		str += " " + ttos(__LINE__) + "\n";
//		str += "Need to set data_ here. I cannot go further until this is implemented\n";
//		throw std::runtime_error(str.c_str());
//	}

	MatrixProductState(IoInputType& io)
		: symmNonconst_(new SymmetryLocalType(io)),symm_(*symmNonconst_)
	{
		// read center and nsites
		io.readline(nsites_,"TotalNumberOfSites=");
		io.readline(center_,"Center=");

		for (size_t i=0;i<nsites_;i++) {
			MpsFactorType f(io,symm_(i));
			data_.push_back(f);
		}
	}

	~MatrixProductState()
	{
		if (symmNonconst_) {
			delete symmNonconst_;
			symmNonconst_=0;
		}
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
	void update(size_t currentSite,const VectorType& v,size_t direction)
	{
		data_[currentSite].updateFromVector(v,direction);
	}

	const MpsFactorType& operator()(size_t site) const
	{
		assert(site<data_.size());
		return data_[site];
	}

	const SymmetryLocalType& symmetry() const { return symm_; }

private:

	// copy ctor:
	MatrixProductState(const MatrixProductState& other);

	// assignment
	MatrixProductState& operator=(const MatrixProductState& other);

	size_t nsites_;
	SymmetryLocalType* symmNonconst_;
	const SymmetryLocalType& symm_;
	size_t center_;
	typename ProgramGlobals::Vector<MpsFactorType>::Type data_;
}; // MatrixProductState

} // namespace Mpspp

/*@}*/
#endif // MATRIX_PRODUCT_STATE_H

