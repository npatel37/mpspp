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

#ifndef CONTRACTED_FACTOR_H
#define CONTRACTED_FACTOR_H

#include "ProgramGlobals.h"

namespace Mpspp {

template<typename MatrixProductOperatorType>
class ContractedFactor {

	typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;
	typedef typename MatrixProductStateType::VectorType VectorType;
	typedef typename MatrixProductStateType::SparseMatrixType SparseMatrixType;
	typedef typename MatrixProductStateType::MpsFactorType  MpsFactorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MatrixProductStateType::ComplexOrRealType ComplexOrRealType;
	typedef ContractedFactor<MatrixProductOperatorType> ThisType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type DenseMatrixType;
	typedef typename MatrixProductStateType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type MatrixType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT, PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:

	typedef typename ProgramGlobals::Vector<SparseMatrixType>::Type DataType;

	ContractedFactor(const MpsFactorType& AorB,
					 const MpoFactorType& h,
					 size_t site,
					 size_t leftOrRight,
					 ThisType* dataPrev,
					 const SymmetryFactorType& symm)
		: data_(h.n_col()),site_(site),leftOrRight_(leftOrRight)
	{
		if (leftOrRight == PART_RIGHT && site>0) data_.resize(h.n_row());

		SparseMatrixType Atranspose;
		transposeConjugate(Atranspose,AorB());
		std::cout<<"ContractedFactor ctor, part="<<partToString();
		std::cout<<" site="<<site<<" size="<<data_.size()<<"\n";
		for (size_t b1=0;b1<data_.size();b1++) {
			if (leftOrRight == PART_RIGHT) {
				initRight2(data_[b1],AorB,b1,h);
			} else {
				initLeft2(data_[b1],AorB,Atranspose,b1,h,dataPrev,symm);
			}
		}
	}

	ContractedFactor(size_t site,
					 size_t leftOrRight)
		: data_(1),site_(site),leftOrRight_(leftOrRight)
	{
		data_[0].makeDiagonal(1,1);
	}

	//! From As (or Bs) and Ws reconstruct *this
	void update(const MpsFactorType& AorB,const MpoFactorType& h,const ThisType& dataPrev,const SymmetryFactorType& symm)
	{
		if (leftOrRight_ == PART_RIGHT) {
			assert(AorB.type()==MpsFactorType::TYPE_B);
			updateRight(AorB,h,dataPrev.data_,symm);
		} else {
			assert(AorB.type()==MpsFactorType::TYPE_A);
			updateLeft(AorB);
		}
	}

	const SparseMatrixType& operator()(size_t b1) const
	{
		assert(b1<data_.size());
		return data_[b1];
	}

	size_t size() const { return data_.size(); }

	template<typename MatrixProductOperatorType_>
	friend std::ostream& operator<<(std::ostream& os,const ContractedFactor<MatrixProductOperatorType_>& contractedFactor);

private:

	void initRight2(SparseMatrixType& m,const MpsFactorType& AorB,size_t b1,const MpoFactorType& h)
	{
		//only when growing:
		contractedFactor0(m,AorB,b1,h);
	}

	void initLeft2old(SparseMatrixType& m,const MpsFactorType& AorB,size_t b,const MpoFactorType& h,ThisType* dataPrev)
	{
		if (site_==0) {
			contractedFactor0(m,AorB,b,h);
			return;
		}

		std::cerr<<"Start initLeft2\n";

		assert(dataPrev!=0);

		const SymmetryFactorType& symm = AorB.symm();
		const SparseMatrixType& A = AorB();
		MatrixType tmp(A.col(),A.col());

		for (size_t b1=0;b1<h.n_row();b1++) {
			const SparseMatrixType& l1 = dataPrev->data_[b1];
			const SparseMatrixType& w = h(b1,b);
			size_t hilbertSize = w.row();
			for (size_t a1=0;a1<l1.row();a1++) {
				for (int k=l1.getRowPtr(a1);k<l1.getRowPtr(a1+1);k++) {
					size_t a1p = l1.getCol(k);
					for (size_t sigma2=0;sigma2<hilbertSize;sigma2++) {
						size_t j = symm.left().pack(a1,sigma2);
						for (int k2=A.getRowPtr(j);k2<A.getRowPtr(j+1);k2++) {
							size_t a2 = A.getCol(k2);
							for (int kp=w.getRowPtr(sigma2);kp<w.getRowPtr(sigma2+1);kp++) {
								size_t sigma2p = w.getCol(kp);
								size_t jp = symm.left().pack(a1p,sigma2p);
								for (int k3=A.getRowPtr(jp);k3<A.getRowPtr(jp+1);k3++) {
									size_t a2p = A.getCol(k3);
									tmp(a2,a2p) += l1.getValue(k)*A.getValue(k2)*w.getValue(kp)*A.getValue(k3);
								} // k3
							} // kp
						} // k2
					} // sigma2
				} // k
			} // a1
		} // b1

		fullMatrixToCrsMatrix(m,tmp);
		std::cerr<<"End initLeft2\n";
	}

	void initLeft2(SparseMatrixType& m,
				   const MpsFactorType& AorB,
				   const SparseMatrixType& Atranspose,
				   size_t b,
				   const MpoFactorType& h,
				   ThisType* dataPrev,
				   const SymmetryFactorType& symm)
	{
		if (site_==0) {
			contractedFactor0(m,AorB,b,h);
			return;
		}

		std::cerr<<"Start initLeft2\n";

		assert(dataPrev!=0);

//		const SymmetryFactorType& symm = AorB.symm();
		const SparseMatrixType& A = AorB();

		m.resize(A.col(),A.col());
		size_t counter=0;
		std::vector<size_t> cols(m.row(),0);
		std::vector<ComplexOrRealType> values(m.row(),0.0);

		for (size_t a2=0;a2<Atranspose.row();a2++) {
			m.setRow(a2,counter);
			for (size_t b1=0;b1<h.n_row();b1++) {
				const SparseMatrixType& l1 = dataPrev->data_[b1];
				const SparseMatrixType& w = h(b1,b);
				if (w.row()==0) continue;
//				size_t hilbertSize = w.row();
				for (int k3=Atranspose.getRowPtr(a2);k3<Atranspose.getRowPtr(a2+1);k3++) {
					PairType a1sigma2 = symm.left().unpack(Atranspose.getCol(k3));
					size_t a1 = a1sigma2.first;
					size_t sigma2 = a1sigma2.second;
					for (int k=l1.getRowPtr(a1);k<l1.getRowPtr(a1+1);k++) {
						size_t a1p = l1.getCol(k);
						for (int kp=w.getRowPtr(sigma2);kp<w.getRowPtr(sigma2+1);kp++) {
							size_t sigma2p = w.getCol(kp);
							size_t j = symm.left().pack(a1p,sigma2p);
							for (int k2=A.getRowPtr(j);k2<A.getRowPtr(j+1);k2++) {
								size_t a2p = A.getCol(k2);
								values[a2p] += l1.getValue(k)*std::conj(Atranspose.getValue(k3))*w.getValue(kp)*A.getValue(k2);
								cols[a2p]=1;
							} // k2
						} // kp
					} // k
				} //k3
			} // b1
			for (size_t i=0;i<cols.size();i++) {
				if (cols[i]==0) continue;
				cols[i]=0;
				m.pushCol(i);
				m.pushValue(values[i]);
				values[i]=0;
				counter++;
			}
		} // a2

		m.setRow(m.row(),counter);
		m.checkValidity();
		std::cerr<<"End initLeft2\n";
	}

	void contractedFactor0(SparseMatrixType& m,const MpsFactorType& AorB,size_t b1,const MpoFactorType& h)
	{
		MatrixType m2(AorB().row(),AorB().col());
		for (size_t a1=0;a1<m2.n_row();a1++) {
			for (size_t a1p=0;a1p<m2.n_row();a1p++) {
				m2(a1,a1p)=contractedFactor0aux(a1,a1p,b1,AorB,h);
			}
		}
		fullMatrixToCrsMatrix(m,m2);
	}

	ComplexOrRealType contractedFactor0aux(size_t a1,size_t a1p,size_t b1,const MpsFactorType& AorB,const MpoFactorType& h) const
	{
		ComplexOrRealType sum = 0;
		MatrixType w;
		size_t row = (leftOrRight_==PART_LEFT) ? 0 : b1;
		size_t col = (leftOrRight_==PART_LEFT) ? b1 : 0;
		if (h.n_row()==1) {
			col=b1;
			row=0;
		}
		crsMatrixToFullMatrix(w,h(row,col));
		for (int k=AorB().getRowPtr(a1);k<AorB().getRowPtr(a1+1);k++) {
			size_t s2 = AorB().getCol(k);
			for (int kp=AorB().getRowPtr(a1p);kp<AorB().getRowPtr(a1p+1);kp++) {
				size_t s2p = AorB().getCol(kp);
				sum += std::conj(AorB().getValue(k)) * w(s2,s2p) * AorB().getValue(kp);
			}
		}
		return sum;
	}

	void updateLeft(const MpsFactorType& A)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to updateLeft(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	void updateRight(const MpsFactorType& B,const MpoFactorType& h,const DataType& dataPrev,const SymmetryFactorType& symm)
	{
		assert(leftOrRight_ == PART_RIGHT);
		assert(B.type()==MpsFactorType::TYPE_B);
		SparseMatrixType Btranspose;
		transposeConjugate(Btranspose,B());
		for (size_t blm2=0;blm2<data_.size();blm2++)
			updateRight(data_[blm2],blm2,B,Btranspose,h,dataPrev,symm);
	}

	void updateRight(SparseMatrixType& m,
					 size_t blm2,
					 const MpsFactorType& B,
					 const SparseMatrixType& Btranspose,
					 const MpoFactorType& h,
					 const DataType& dataPrev,
					 const SymmetryFactorType& symm)
	{
		size_t someSize = B().row();
		m.resize(someSize,someSize);

		size_t counter = 0;
		std::vector<size_t> cols(m.row(),0);
		std::vector<ComplexOrRealType> values(m.row(),0.0);

		for (size_t alm2=0;alm2<m.row();alm2++) {
			m.setRow(alm2,counter);
			updateRight(values,cols,alm2,blm2,B,Btranspose,h,dataPrev,symm);
			for (size_t i=0;i<cols.size();i++) {
				if (cols[i]==0) continue;
				cols[i]=0;
				m.pushCol(i);
				m.pushValue(values[i]);
				values[i]=0;
				counter++;
			}
		}
		m.setRow(m.row(),counter);
	}

	void updateRight(std::vector<ComplexOrRealType>& values,
					 std::vector<size_t>& cols,
					 size_t alm2,
					 size_t blm2,
					 const MpsFactorType& B,
					 const SparseMatrixType& Btranspose,
					 const MpoFactorType& h,
					 const DataType& dataPrev,
					 const SymmetryFactorType& symm)
	{
//		const SymmetryFactorType& symm = B.symm();
		const SparseMatrixType& Bmatrix = B();

		for (int kb=Bmatrix.getRowPtr(alm2);kb<Bmatrix.getRowPtr(alm2+1);kb++) {
			PairType sigmalm1alm1 = symm.right().unpack(Bmatrix.getCol(kb));
			size_t sigmalm1 = sigmalm1alm1.first;
			size_t alm1 = sigmalm1alm1.second;
			for (size_t blm1=0;blm1<dataPrev.size();blm1++) {
				const SparseMatrixType& w = h(blm2,blm1);
				if (w.row()==0) continue;
				for (int k=w.getRowPtr(sigmalm1);k<w.getRowPtr(sigmalm1+1);k++) {
					size_t sigmaplm1 = w.getCol(k);
					for (int kprev=dataPrev[blm1].getRowPtr(alm1);kprev<dataPrev[blm1].getRowPtr(alm1+1);kprev++) {
						size_t aplm1 = dataPrev[blm1].getCol(kprev);
						size_t i = symm.right().pack(sigmaplm1,aplm1);
						for (int kb2=Btranspose.getRowPtr(i);kb2<Btranspose.getRowPtr(i+1);kb2++) {
							size_t aplm2 = Btranspose.getCol(kb2);
							values[aplm2] += std::conj(Bmatrix.getValue(kb))*w.getValue(k)*Btranspose.getValue(kb2) * dataPrev[blm1].getValue(kprev);
							cols[aplm2]=1;
						} // kb2
					} // kprev
				} // k
			} // blm1
		} // kb
	}

	std::string partToString() const
	{
		return (leftOrRight_==PART_LEFT) ? "left" : "right";
	}

	DataType data_;
	size_t site_;
	size_t leftOrRight_;

}; // ContractedFactor

template<typename MatrixProductOperatorType>
std::ostream& operator<<(std::ostream& os,const ContractedFactor<MatrixProductOperatorType>& contractedFactor)
{
	os<<"SparseMatrices= "<<contractedFactor.data_.size()<<"\n";
	for (size_t i=0;i<contractedFactor.data_.size();i++) {
		os<<contractedFactor.data_[i];
	}
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_FACTOR_H

