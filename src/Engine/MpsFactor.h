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

	typedef std::pair<SizeType,SizeType> PairType;

public:

	enum {TYPE_A,TYPE_B};

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef VectorWithOffset<ComplexOrRealType> VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename SymmetryFactorType::IoInputType IoInputType;
	typedef typename SymmetryFactorType::SymmetryComponentType SymmetryComponentType;
	typedef typename SymmetryComponentType::VectorIntegerType VectorIntegerType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::RandomForTests<RealType> RandomNumberGeneratorType;

	MpsFactor(SizeType aOrB)
	    : rng_(0),aOrB_(aOrB)
	{}

	void setRandom(SizeType site,SizeType n)
	{
		data_.resize(n,n);
		data_.makeDiagonal(n,1.0);

		//assert(isNormalized(data_));
	}

	template<typename SomeTruncationType>
	void move(SomeTruncationType& truncation,
	          const VectorType& v,
	          SizeType symmetrySector,
	          const SymmetryFactorType& symm)
	{

		SizeType row = symm.left().size();
		SizeType col = symm.right().size();
		if (aOrB_==TYPE_B) std::swap(row,col);

		MatrixType m(row,col);

		SizeType offset = symm.super().partitionOffset(symmetrySector);
		SizeType total = symm.super().partitionSize(symmetrySector);
		for (SizeType i=0;i<total;i++) {
			PairType ab = symm.super().unpack(i+offset);
			if (aOrB_==TYPE_A) {
				m(ab.first,ab.second) = v[i];
			} else {
				m(ab.second,ab.first) = v[i];
			}
		}

		moveFromVector(m,truncation,symm,symmetrySector);
	}

	template<typename SomeTruncationType>
	void truncate(SizeType cutoff,const SomeTruncationType& trunc)
	{
		trunc.matrixRow(data_,cutoff);
	}

	PsimagLite::String typeToString() const
	{
		return (aOrB_==TYPE_A) ? "A" : "B";
	}

	const SparseMatrixType& operator()() const { return data_; }

	const SizeType type() const { return aOrB_; }

	template<typename ComplexOrRealType2,typename SymmetryLocalType2>
	friend std::ostream& operator<<(std::ostream&,
	                                const MpsFactor<ComplexOrRealType2,SymmetryLocalType2>&);

private:

	template<typename SomeTruncationType>
	void moveFromVector(const MatrixType& m,
	                    SomeTruncationType& truncation,
	                    const SymmetryFactorType& symm,
	                    SizeType)
	{
		const SymmetryComponentType& summed = (aOrB_==TYPE_A) ? symm.left() : symm.right();
		const SymmetryComponentType& nonSummed = (aOrB_==TYPE_A) ? symm.right() : symm.left();

		MatrixType finalU(summed.size(),summed.size());

		truncation.setSize(summed.size());
		for (SizeType i=0;i<summed.partitions()-1;i++) {
			SizeType istart = summed.partitionOffset(i);
			SizeType itotal = summed.partitionSize(i);
			for (SizeType j=0;j<nonSummed.partitions()-1;j++) {
				SizeType jstart = nonSummed.partitionOffset(j);
				SizeType jtotal = nonSummed.partitionSize(j);

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
		// debuggin only
		assert(aOrB_ == TYPE_A);
		SizeType n = data_.row();
		data_.makeDiagonal(n,1.0);
	}

	void setThisSector(MatrixType& u,
	                   SizeType istart,
	                   SizeType itotal,
	                   SizeType jstart,
	                   SizeType jtotal,
	                   const MatrixType& m) const
	{
		for (SizeType i=0;i<itotal;i++) {
			for (SizeType j=0;j<jtotal;j++) {
				u(i,j) = m(i+istart,j+jstart);
			}
		}
	}

	void setFinalU(MatrixType& finalU,
	               SizeType istart,
	               SizeType itotal,
	               SizeType jstart,
	               SizeType jtotal,
	               const MatrixType& u) const
	{
		// std::cout<<"setFinalU from "<<istart<<" to "<<(istart+itotal-1)<<"\n";
		SizeType n = std::min(itotal,u.n_col());
		for (SizeType i=0;i<itotal;i++) {
			for (SizeType j=0;j<n;j++) {
				if (j+istart>=finalU.n_col()) continue;
				finalU(i+istart,j+istart) = u(i,j);
			}
		}
	}

	template<typename SomeTruncationType>
	void setFinalS(SomeTruncationType& truncation,
	               SizeType jstart,
	               SizeType jtotal,
	               const VectorRealType& s) const
	{
		//	std::cout<<"setFinalS from "<<jstart<<" to "<<(jstart+jtotal-1)<<"\n";
		SizeType n = std::min(jtotal,static_cast<SizeType>(s.size()));
		for (SizeType j=0;j<n;j++) {
			//			assert(j+jstart<finalS.size());
			//			if (fabs(s[j])<1e-6) continue;
			if (j+jstart>=truncation.size()) continue;
			truncation(j+jstart) = s[j];
			//			std::cout<<"s["<<j<<"]="<<s[j]<<" ";
		}
		//		std::cout<<"\n";
	}

	void setFinalVt(MatrixType& finalVt,
	                SizeType istart,
	                SizeType itotal,
	                SizeType jstart,
	                SizeType jtotal,
	                const MatrixType& vt) const
	{
		std::cout<<"setFinalU from "<<istart<<" to "<<(istart+itotal-1)<<"\n";
		for (SizeType i=0;i<jtotal;i++) {
			for (SizeType j=0;j<jtotal;j++) {
				finalVt(i+jstart,j+jstart) = vt(i,j);
			}
		}
	}

	template<typename SomeTruncationType>
	bool isCorrectSvd(const MatrixType& mat,
	                  const MatrixType& u,
	                  SomeTruncationType& truncation,
	                  const MatrixType& vt) const
	{
		MatrixType m(u.n_row(),vt.n_col());
		truncation.recoverSvd(m,u,vt);
		std::cout<<"recovering matrix from svd:\n";
		std::cout<<m;
		for (SizeType i=0;i<m.n_row();i++) {
			for (SizeType j=0;j<m.n_col();j++)
				if (fabs(m(i,j)-mat(i,j))>1e-6) return false;
		}
		return true;
	}

	bool respectsSymmetry(const MatrixType& m,
	                      const SymmetryComponentType& summed) const
	{
		assert(m.n_row()<=summed.size());
		assert(m.n_col()<=summed.size());
		for (SizeType i=0;i<m.n_row();i++) {
			SizeType qi = summed.qn(i);
			for (SizeType j=0;j<m.n_col();j++) {
				SizeType qj = summed.qn(j);
				if (qi==qj) continue;
				if (fabs(m(i,j))>1e-6) return false;
			}
		}
		return true;
	}

	bool isNormalized(const MatrixType& m) const
	{
		//		assert(m.n_row()==m.n_col());
		for (SizeType i=0;i<m.n_col();i++) {
			for (SizeType j=0;j<m.n_col();j++) {
				ComplexOrRealType sum = 0;
				for (SizeType k=0;k<m.n_row();k++) {
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
	SizeType aOrB_;
}; // MpsFactor

template<typename ComplexOrRealType,typename SymmetryLocalType>
std::ostream& operator<<(std::ostream& os,
                         const MpsFactor<ComplexOrRealType,SymmetryLocalType>& mps)
{
	os<<"type= "<<mps.typeToString();
	os<<" rows= "<<mps.data_.row()<<" cols="<<mps.data_.col()<<"\n";
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // MPS_FACTOR_TYPE_H

