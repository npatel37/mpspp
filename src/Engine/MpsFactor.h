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

#ifndef MPS_FACTOR_TYPE_H
#define MPS_FACTOR_TYPE_H

#include "VectorWithOffset.h"
#include "ProgramGlobals.h"
#include "RandomForTests.h"

namespace Mpspp {

template<typename ComplexOrRealType,typename SymmetryFactorType>
class MpsFactor {

	typedef std::pair<size_t,size_t> PairType;

public:

	enum {TYPE_A,TYPE_B};

	typedef typename ProgramGlobals::Real<ComplexOrRealType>::Type RealType;
	typedef VectorWithOffset<ComplexOrRealType> VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename ProgramGlobals::CrsMatrix<ComplexOrRealType>::Type SparseMatrixType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type MatrixType;
	typedef typename SymmetryFactorType::IoInputType IoInputType;
	typedef PsimagLite::RandomForTests<RealType> RandomNumberGeneratorType;

	MpsFactor(const SymmetryFactorType& symm,size_t site,size_t aOrB)
		: symm_(symm),site_(site),rng_(0),aOrB_(aOrB)
	{}

//	MpsFactor(IoInputType& io,const SymmetryFactorType& symm,size_t site)
//	: symm_(symm),rng_(0)
//	{
////		if (site==0 || site+1 ==symm_.super().block().size()) {
////			std::string str(__FILE__);
////			str += " " + ttos(__LINE__) + "\n";
////			str += "Needs corner cases. I cannot go further until this is implemented\n";
////			std::cerr<<str;
////			return;
////		}
//		VectorWithOffsetType mpsFactorM;
//		mpsFactorM.load(io,"MpsFactorM");

//		MatrixType mpsFactor;
//		findMpsFactor1(mpsFactor,mpsFactorM);

//		std::vector<RealType> s;
//		svd(mpsFactor,s);

//		fullMatrixToCrsMatrix(data_,mpsFactor);
//	}

	MpsFactor& operator=(const MpsFactor& other)
	{
		this->data_ = other.data_;
		this->site_ = other.site_;
//		this->rng_ = other.rng_;
		this->aOrB_ = other.aOrB_;
		//this->symm_ = other.symm_;
		return *this;
	}

	void setRandom(size_t site)
	{
		MatrixType m(symm_.right().size(),symm_.right().size());
		for (size_t i=0;i<m.n_row();i++)
			for (size_t j=0;j<m.n_col();j++)
				m(i,j) = rng_();
		fullMatrixToCrsMatrix(data_,m);
	}

	void updateFromVector(const VectorType& v,size_t symmetrySector)
	{
		MatrixType m(symm_.left().size(),symm_.right().size());


		size_t offset = symm_.super().partitionOffset(symmetrySector);
		size_t total = symm_.super().partitionSize(symmetrySector);
		for (size_t i=0;i<total;i++) {
			PairType ab = symm_.super().unpack(i+offset);
			if (aOrB_==TYPE_A) {
				m(ab.first,ab.second) = v[i];
			} else {
				m(ab.second,ab.first) = v[i];
			}
		}

		std::vector<RealType> s;
		svd(m,s);
		fullMatrixToCrsMatrix(data_,m);
	}

	const SparseMatrixType& operator()() const { return data_; }

	const SymmetryFactorType& symm() const { return symm_; }

	const size_t type() const { return aOrB_; }

private:

	const SymmetryFactorType& symm_;
	size_t site_;
	RandomNumberGeneratorType rng_;
	SparseMatrixType data_;
	size_t aOrB_;
}; // MpsFactor

} // namespace Mpspp

/*@}*/
#endif // MPS_FACTOR_TYPE_H

