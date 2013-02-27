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

#ifndef MPS_LOCAL_H
#define MPS_LOCAL_H

#include "ProgramGlobals.h"
#include "MpsFactor.h"

namespace Mpspp {

template<typename ComplexOrRealType_,typename SymmetryLocalType_>
class MpsLocal {

	// FIXME: IDEA: PULL SYMMETRY OUT, PASS THROUGH FUNCTIONS

public:

	typedef SymmetryLocalType_ SymmetryLocalType;
	typedef typename SymmetryLocalType::IoInputType IoInputType;
	typedef ComplexOrRealType_ ComplexOrRealType;
	typedef typename ProgramGlobals::Real<ComplexOrRealType>::Type RealType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef MpsFactor<ComplexOrRealType,SymmetryFactorType> MpsFactorType;
	typedef typename MpsFactorType::VectorType VectorType;
	typedef typename MpsFactorType::SparseMatrixType SparseMatrixType;
	typedef typename MpsFactorType::MatrixType MatrixType;
	typedef typename MpsFactorType::VectorRealType VectorRealType;
	typedef typename MpsFactorType::VectorIntegerType VectorIntegerType;

	MpsLocal(size_t nsites)
	: nsites_(nsites),center_(0)
	{}

//	MpsLocal(IoInputType& io)
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

	~MpsLocal()
	{
		for (size_t i=0;i<A_.size();i++)
			if (A_[i]) delete A_[i];
		for (size_t i=0;i<B_.size();i++)
			if (B_[i]) delete B_[i];
	}

	void grow(size_t currentSite,const SymmetryLocalType& symm,size_t nk)
	{
		center_ = currentSite;
//		if (currentSite==0) {
//			MpsFactorType* mpsFactor2 = new MpsFactorType(MpsFactorType::TYPE_B);
//			mpsFactor2->setRandom(currentSite,nk);
//			B_.push_back(mpsFactor2);
//		}

		MpsFactorType* mpsFactor = new MpsFactorType(MpsFactorType::TYPE_B);
		size_t n = symm(currentSite+1).right().size();
		mpsFactor->setRandom(currentSite,n);
		B_.push_back(mpsFactor);


//		if (currentSite==0) {
//			MpsFactorType* mpsFactor3 = new MpsFactorType(MpsFactorType::TYPE_A);
//			mpsFactor3->setRandom(currentSite,nk);
//			A_.push_back(mpsFactor3);
//		}
		MpsFactorType* mpsFactor2 = new MpsFactorType(MpsFactorType::TYPE_A);
		n =  symm(currentSite+1).left().size();
		mpsFactor2->setRandom(currentSite,n);
		A_.push_back(mpsFactor2);
	}

//	size_t center() const
//	{
//		return center_;
//	}

	//! Returns the number of sites
	size_t sites() const
	{
		return nsites_;
	}

	//! tmpVec[i] --> M^\sigma2 _ {a1,a2}
	template<typename SomeTruncationType>
	void move(SomeTruncationType& truncation,
			  size_t currentSite,
			  const VectorType& v,
			  size_t direction,
			  size_t symmetrySector,
			  const SymmetryFactorType& symm)
	{
		center_ = currentSite;
		if (direction==ProgramGlobals::TO_THE_RIGHT) {
			if (currentSite==A_.size()) {
				MpsFactorType* mpsFactor = new MpsFactorType(MpsFactorType::TYPE_A);
				size_t n = symm.left().size();
				mpsFactor->setRandom(currentSite,n);
				A_.push_back(mpsFactor);
			}
			assert(currentSite<A_.size());
			A_[currentSite]->move(truncation,v,symmetrySector,symm);
			std::cout<<"moved A["<<currentSite<<"].row= ";
			std::cout<<A_[currentSite]->operator()().row()<<"\n";
		} else {
			size_t nsites = symm.super().block().size();
			size_t siteToSet = nsites - 1 - currentSite;
			assert(siteToSet<nsites);

			if (siteToSet==B_.size()) {
				MpsFactorType* mpsFactor = new MpsFactorType(MpsFactorType::TYPE_B);
				size_t n = symm.right().size();
				mpsFactor->setRandom(currentSite,n);
				B_.push_back(mpsFactor);
			}
			B_[siteToSet]->move(truncation,v,symmetrySector,symm);
			std::cout<<"moved B["<<(siteToSet)<<"].row= ";
			std::cout<<B_[siteToSet]->operator()().row()<<"\n";
		}
	}

	template<typename SomeTruncationType>
	void truncate(size_t site,size_t part,size_t cutoff,size_t nsites,const SomeTruncationType& trunc)
	{
		if (part==ProgramGlobals::PART_LEFT) {
			A_[site]->truncate(cutoff,trunc);
		} else {
			size_t siteToSet = nsites-site-1;
			B_[siteToSet]->truncate(cutoff,trunc);
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


	typename ProgramGlobals::Real<ComplexOrRealType>::Type norm(size_t type,const SymmetryLocalType& symm) const
	{
		if (type==MpsFactorType::TYPE_B) {
			return normB(symm);
		}
		return normA(symm);
	}

	template<typename ComplexOrRealType2,typename SymmetryLocalType2>
	friend std::ostream& operator<<(std::ostream& os,const MpsLocal<ComplexOrRealType2,SymmetryLocalType2>& mps);

private:

	typename ProgramGlobals::Real<ComplexOrRealType>::Type normB(const SymmetryLocalType& symm) const
	{
		MatrixType tmpOld;

		for (size_t i=0;i<B_.size();i++) {
			size_t center = B_.size()-1-i;
			MatrixType tmpNew;
			computeIntermediate(tmpNew,tmpOld,center,symm);
			tmpOld = tmpNew;
		}

		ComplexOrRealType sum = 0;
		for (size_t i=0;i<tmpOld.n_row();i++) {
			for (size_t j=0;j<tmpOld.n_col();j++) {
				sum += tmpOld(i,j);
			}
		}
		return std::real(sum);
	}

	void computeIntermediate(MatrixType& matrixNew,const MatrixType& matrixOld,size_t center,const SymmetryLocalType& symmLocal) const
	{
		const SparseMatrixType& Bmatrix = B_[center]->operator()();
		SparseMatrixType Btranspose;
		transposeConjugate(Btranspose,Bmatrix);

		matrixNew.resize(Bmatrix.row(),Bmatrix.row());
		matrixNew.setTo(0.0);

		if (matrixOld.n_row()==0) {
			assert(center+1==B_.size());
			for (size_t x=0;x<Bmatrix.row();x++) {
				for (int k=Bmatrix.getRowPtr(x);k<Bmatrix.getRowPtr(x+1);k++) {
					size_t sigma = Bmatrix.getCol(k);
					for (int k2=Btranspose.getRowPtr(sigma);k2<Btranspose.getRowPtr(sigma+1);k2++) {
						size_t y = Btranspose.getCol(k2);
						matrixNew(x,y) += Bmatrix.getValue(k) * std::conj(Btranspose.getValue(k2));
					}
				}
			}
			return;
		}

		const SymmetryFactorType& symm = symmLocal(center);
		for (size_t anm2=0;anm2<Bmatrix.row();anm2++) {
			for (int k=Bmatrix.getRowPtr(anm2);k<Bmatrix.getRowPtr(anm2+1);k++) {
				PairType p = symm.right().unpack(Bmatrix.getCol(k));
				size_t sigmanm2=p.first;
				size_t anm1=p.second;
				for (size_t apnm1=0;apnm1<matrixOld.n_row();apnm1++) {
					size_t x = symm.right().pack(sigmanm2,apnm1);
					for (int k2=Btranspose.getRowPtr(x);k2<Btranspose.getRowPtr(x+1);k2++) {
						size_t apnm2 = Btranspose.getCol(k2);
						matrixNew(anm2,apnm2) += matrixOld(anm1,apnm1) * Bmatrix.getValue(k) * Btranspose.getValue(k2);
					} // k2
				} // apnm1
			} // k
		} // anm2
	}

	typename ProgramGlobals::Real<ComplexOrRealType>::Type normA(const SymmetryLocalType& symm) const
	{
		MatrixType tmpOld;

		for (size_t i=0;i<A_.size();i++) {
			size_t center = i;
			MatrixType tmpNew;
			computeIntermediateA(tmpNew,tmpOld,center,symm);
			tmpOld = tmpNew;
		}

		ComplexOrRealType sum = 0;
		for (size_t i=0;i<tmpOld.n_row();i++) {
			for (size_t j=0;j<tmpOld.n_col();j++) {
				sum += tmpOld(i,j);
			}
		}
		return std::real(sum);
	}

	void computeIntermediateA(MatrixType& matrixNew,const MatrixType& matrixOld,size_t center,const SymmetryLocalType& symmLocal) const
	{
		const SymmetryFactorType& symm = symmLocal(center+1);
		size_t leftSize = symm.left().split();
		size_t leftSize2 = symm.left().size();

		SparseMatrixType Atranspose;
		transposeConjugate(Atranspose,A_[center]->operator()());
		matrixNew.resize(leftSize*leftSize,leftSize2*leftSize2);
		matrixNew.setTo(0.0);
		assert(Atranspose.row()==symm.left().size());
		for (size_t a1=0;a1<Atranspose.row();a1++) {
			for (int k1=Atranspose.getRowPtr(a1);k1<Atranspose.getRowPtr(a1+1);k1++) {
				PairType a0sigma1=symm.left().unpack(Atranspose.getCol(k1));
				size_t a0 = a0sigma1.first;
				size_t sigma1= a0sigma1.second;
				for (size_t a1p=0;a1p<Atranspose.row();a1p++) {
					for (int k2=Atranspose.getRowPtr(a1p);k2<Atranspose.getRowPtr(a1p+1);k2++) {

						PairType a0sigma1p=symm.left().unpack(Atranspose.getCol(k2));
						size_t a0p = a0sigma1p.first;
						size_t sigma1p= a0sigma1p.second;
						if (sigma1!=sigma1p) continue;
						matrixNew(a0+a0p*leftSize,a1+a1p*leftSize2) += Atranspose.getValue(k1) * std::conj(Atranspose.getValue(k2));
					} // k2
				} // a1p
			} // k1
		} // a1

		if (matrixOld.n_row()==0) return;
		MatrixType m = matrixOld * matrixNew;
		matrixNew = m;
	}

	// copy ctor:
	MpsLocal(const MpsLocal& other);

	// assignment
	MpsLocal& operator=(const MpsLocal& other);

//	const SymmetryLocalType& symm_;
	size_t nsites_;
	size_t center_;
	typename ProgramGlobals::Vector<MpsFactorType*>::Type B_;
	typename ProgramGlobals::Vector<MpsFactorType*>::Type A_;
}; // MpsLocal

template<typename ComplexOrRealType,typename SymmetryLocalType>
std::ostream& operator<<(std::ostream& os,const MpsLocal<ComplexOrRealType,SymmetryLocalType>& mps)
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
#endif // MPS_LOCAL_H

