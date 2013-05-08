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
#include "Sort.h"

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
	typedef typename SymmetryFactorType::SymmetryComponentType SymmetryComponentType;
	typedef typename SymmetryComponentType::VectorIntegerType VectorIntegerType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::RandomForTests<RealType> RandomNumberGeneratorType;

	MpsFactor(size_t aOrB)
		: rng_(0),aOrB_(aOrB)
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

//		VectorRealType s;
//		svd(mpsFactor,s);

//		fullMatrixToCrsMatrix(data_,mpsFactor);
//	}

	void setRandom(size_t site,size_t n)
	{
//		MatrixType m(n,n);
//		for (size_t i=0;i<n;i++)
//			m(i,i) = (rng_()<0.5) ? -1 : 1;

//		fullMatrixToCrsMatrix(data_,m);
		data_.resize(n,n);
		data_.makeDiagonal(n,1.0);

		assert(isNormalized(data_));
	}

	template<typename SomeTruncationType>
	void move(SomeTruncationType& truncation,const VectorType& v,size_t symmetrySector,const SymmetryFactorType& symm)
	{

		size_t row = symm.left().size();
		size_t col = symm.right().size();
		if (aOrB_==TYPE_B) std::swap(row,col);

		MatrixType m(row,col);

//		size_t total = symm.super().size();
		size_t offset = symm.super().partitionOffset(symmetrySector);
		size_t total = symm.super().partitionSize(symmetrySector);
		for (size_t i=0;i<total;i++) {
			PairType ab = symm.super().unpack(i+offset);
			if (aOrB_==TYPE_A) {
				m(ab.first,ab.second) = v[i];
			} else {
				m(ab.second,ab.first) = v[i];
			}
		}
//		std::cout<<"full matrix prior to svd is\n";
//		std::cout<<m;

		moveFromVector(m,truncation,symm,symmetrySector);
	}

	template<typename SomeTruncationType>
	void truncate(size_t cutoff,const SomeTruncationType& trunc)
	{
		trunc.matrixRow(data_,cutoff);
	}

//	template<typename SomeNumericType>
//	void divideBy(const SomeNumericType& value)
//	{
//		assert(value!=0);
//		data_ *= (1.0/value);
//	}

	std::string typeToString() const
	{
		return (aOrB_==TYPE_A) ? "A" : "B";
	}

	const SparseMatrixType& operator()() const { return data_; }

	const size_t type() const { return aOrB_; }

	template<typename ComplexOrRealType2,typename SymmetryLocalType2>
	friend std::ostream& operator<<(std::ostream& os,const MpsFactor<ComplexOrRealType2,SymmetryLocalType2>& mps);

private:

	template<typename SomeTruncationType>
	void moveFromVector(const MatrixType& m,SomeTruncationType& truncation,const SymmetryFactorType& symm,size_t symmetrySector)
	{
		const SymmetryComponentType& summed = (aOrB_==TYPE_A) ? symm.left() : symm.right();
		const SymmetryComponentType& nonSummed = (aOrB_==TYPE_A) ? symm.right() : symm.left();

//		size_t max = std::max(summed.size(),nonSummed.size());
		MatrixType finalU(summed.size(),summed.size());

		truncation.setSize(summed.size());
		for (size_t i=0;i<summed.partitions()-1;i++) {
			size_t istart = summed.partitionOffset(i);
			size_t itotal = summed.partitionSize(i);
			for (size_t j=0;j<nonSummed.partitions()-1;j++) {
				size_t jstart = nonSummed.partitionOffset(j);
				size_t jtotal = nonSummed.partitionSize(j);

				MatrixType u(itotal,jtotal);
				setThisSector(u,istart,itotal,jstart,jtotal,m);

				VectorRealType s;
				MatrixType vt;
				svd('A',u,s,vt);

				setFinalU(finalU,istart,itotal,jstart,jtotal,u);
				setFinalS(truncation,istart,itotal,s);
			}
		}

//		VectorRealType finalS(m.n_col());
//		svd('A',finalU,finalS,finalVt);
//		truncation.set(finalS);

		MatrixType mtranspose;
		if (aOrB_==TYPE_B)
			transposeConjugate(mtranspose,finalU);

//		std::cout<<"new AorB=\n";
//		std::cout<<finalU;
//		truncation.print(std::cout);
//		std::cout<<"final vt\n";
//		std::cout<<finalVt;
		assert(isNormalized(finalU));
		assert(respectsSymmetry(finalU,summed));
//		assert(isCorrectSvd(m,finalU,truncation,finalVt));
		fullMatrixToCrsMatrix(data_,(aOrB_==TYPE_A) ? finalU : mtranspose);
	}

	void setThisSector(MatrixType& u,
					   size_t istart,
					   size_t itotal,
					   size_t jstart,
					   size_t jtotal,
					   const MatrixType& m) const
	{
		for (size_t i=0;i<itotal;i++) {
			for (size_t j=0;j<jtotal;j++) {
				u(i,j) = m(i+istart,j+jstart);
			}
		}
	}

	void setFinalU(MatrixType& finalU,
				   size_t istart,
				   size_t itotal,
				   size_t jstart,
				   size_t jtotal,
				   const MatrixType& u) const
	{
//		std::cout<<"setFinalU from "<<istart<<" to "<<(istart+itotal-1)<<"\n";
		size_t n = std::min(itotal,u.n_col());
		for (size_t i=0;i<itotal;i++) {
			for (size_t j=0;j<n;j++) {
				if (j+istart>=finalU.n_col()) continue;
				finalU(i+istart,j+istart) = u(i,j);
			}
		}
	}

	template<typename SomeTruncationType>
	void setFinalS(SomeTruncationType& truncation,
				   size_t jstart,
				   size_t jtotal,
				   const VectorRealType& s) const
	{
//		std::cout<<"setFinalS from "<<jstart<<" to "<<(jstart+jtotal-1)<<"\n";
		size_t n = std::min(jtotal,s.size());
		for (size_t j=0;j<n;j++) {
//			assert(j+jstart<finalS.size());
//			if (fabs(s[j])<1e-6) continue;
			if (j+jstart>=truncation.size()) continue;
			truncation(j+jstart) = s[j];
//			std::cout<<"s["<<j<<"]="<<s[j]<<" ";
		}
//		std::cout<<"\n";
	}

	void setFinalVt(MatrixType& finalVt,
				   size_t istart,
				   size_t itotal,
				   size_t jstart,
				   size_t jtotal,
				   const MatrixType& vt) const
	{
		std::cout<<"setFinalU from "<<istart<<" to "<<(istart+itotal-1)<<"\n";
		for (size_t i=0;i<jtotal;i++) {
			for (size_t j=0;j<jtotal;j++) {
				finalVt(i+jstart,j+jstart) = vt(i,j);
			}
		}
	}

	template<typename SomeTruncationType>
	bool isCorrectSvd(const MatrixType& mat,const MatrixType& u,SomeTruncationType& truncation,const MatrixType& vt) const
	{
		MatrixType m(u.n_row(),vt.n_col());
		truncation.recoverSvd(m,u,vt);
		std::cout<<"recovering matrix from svd:\n";
		std::cout<<m;
		for (size_t i=0;i<m.n_row();i++) {
			for (size_t j=0;j<m.n_col();j++)
				if (fabs(m(i,j)-mat(i,j))>1e-6) return false;
		}
		return true;
	}

	bool respectsSymmetry(const MatrixType& m,const SymmetryComponentType& summed) const
	{
		assert(m.n_row()<=summed.size());
		assert(m.n_col()<=summed.size());
		for (size_t i=0;i<m.n_row();i++) {
			size_t qi = summed.qn(i);
			for (size_t j=0;j<m.n_col();j++) {
				size_t qj = summed.qn(j);
				if (qi==qj) continue;
				if (fabs(m(i,j))>1e-6) return false;
			}
		}
		return true;
	}

	bool isNormalized(const MatrixType& m) const
	{
//		assert(m.n_row()==m.n_col());
		for (size_t i=0;i<m.n_col();i++) {
			for (size_t j=0;j<m.n_col();j++) {
				ComplexOrRealType sum = 0;
				for (size_t k=0;k<m.n_row();k++) {
					sum += m(k,i) * std::conj(m(k,j));
				}
				if (i==j && fabs(sum-1.0)>1e-5)
					return false;
				if (i!=j && fabs(sum)>1e-5)
					return false;
			}
		}
		return true;
	}

	RandomNumberGeneratorType rng_;
	SparseMatrixType data_;
	size_t aOrB_;
}; // MpsFactor

template<typename ComplexOrRealType,typename SymmetryLocalType>
std::ostream& operator<<(std::ostream& os,const MpsFactor<ComplexOrRealType,SymmetryLocalType>& mps)
{
	os<<"type= "<<mps.typeToString();
	os<<" rows= "<<mps.data_.row()<<" cols="<<mps.data_.col()<<"\n";
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // MPS_FACTOR_TYPE_H

