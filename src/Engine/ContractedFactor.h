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

	typedef typename MatrixProductOperatorType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::VectorType VectorType;
	typedef typename MpsLocalType::SparseMatrixType SparseMatrixType;
	typedef typename MpsLocalType::MpsFactorType  MpsFactorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MpsLocalType::ComplexOrRealType ComplexOrRealType;
	typedef ContractedFactor<MatrixProductOperatorType> ThisType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type DenseMatrixType;
	typedef typename MpsLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type MatrixType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT, PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:

	typedef typename ProgramGlobals::Vector<SparseMatrixType>::Type DataType;

	void build(const MpsFactorType& AorB,const MpoFactorType& h,const ThisType& prev,const SymmetryFactorType& symm)
	{
		if (leftOrRight_==PART_RIGHT) {
			assert(AorB.type()==MpsFactorType::TYPE_B);

			data_.resize(h.n_row());
			moveRight(AorB,h,prev.data_,symm);
		} else {
			assert(AorB.type()==MpsFactorType::TYPE_A);

			data_.resize(h.n_col());
			moveLeft(AorB,h,prev.data_,symm);
		}
	}

//	void build(const MpoFactorType& h)
//	{
//		assert(leftOrRight_==PART_LEFT);
//		data_.resize(h.n_row());
//	}

	ContractedFactor(size_t leftOrRight)
		: data_(1),leftOrRight_(leftOrRight)
	{
		data_[0].makeDiagonal(1,1);
	}

	//! From As (or Bs) and Ws reconstruct *this
	void move(const MpsFactorType& AorB,const MpoFactorType& h,const ThisType& dataPrev,const SymmetryFactorType& symm)
	{
		if (leftOrRight_ == PART_RIGHT) {
			assert(AorB.type()==MpsFactorType::TYPE_B);
			moveRight(AorB,h,dataPrev.data_,symm);
		} else {
			assert(AorB.type()==MpsFactorType::TYPE_A);
			moveLeft(AorB,h,dataPrev.data_,symm);
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

	void moveLeft(SparseMatrixType& m,
				   const MpsFactorType& AorB,
				   const SparseMatrixType& Atranspose,
				   size_t b,
				   const MpoFactorType& h,
				   const DataType& dataPrev,
				   const SymmetryFactorType& symm)
	{
		std::cerr<<"Start initLeft2\n";

		const SparseMatrixType& A = AorB();

		m.resize(A.col(),A.col());
		size_t counter=0;
		std::vector<size_t> cols(m.row(),0);
		std::vector<ComplexOrRealType> values(m.row(),0.0);
		assert(symm.left().split()==0 || symm.left().split()==dataPrev[0].row());
		assert(symm.left().size()==A.row());
		for (size_t a2=0;a2<Atranspose.row();a2++) {
			m.setRow(a2,counter);
			for (int k3=Atranspose.getRowPtr(a2);k3<Atranspose.getRowPtr(a2+1);k3++) {
				PairType a1sigma2 = symm.left().unpack(Atranspose.getCol(k3));
				size_t a1 = a1sigma2.first;
				size_t sigma2 = a1sigma2.second;

				for (size_t b1=0;b1<h.n_row();b1++) {
					const SparseMatrixType& w = h(b1,b);
					if (w.row()==0) continue;
					const SparseMatrixType& l1 = dataPrev[b1];

					for (int k=l1.getRowPtr(a1);k<l1.getRowPtr(a1+1);k++) {
						size_t a1p = l1.getCol(k);
						ComplexOrRealType tmp= l1.getValue(k)*std::conj(Atranspose.getValue(k3));
						for (int kp=w.getRowPtr(sigma2);kp<w.getRowPtr(sigma2+1);kp++) {
							size_t sigma2p = w.getCol(kp);
							size_t j = symm.left().pack(a1p,sigma2p);
							for (int k2=A.getRowPtr(j);k2<A.getRowPtr(j+1);k2++) {
								size_t a2p = A.getCol(k2);
								values[a2p] += tmp*w.getValue(kp)*A.getValue(k2);
								cols[a2p]=1;
							} // k2
						} // kp
					} // k
				} //b1
			} // k3
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

	void moveLeft(const MpsFactorType& A,const MpoFactorType& h,const DataType& dataPrev,const SymmetryFactorType& symm)
	{
		assert(leftOrRight_ == PART_LEFT);
		assert(A.type()==MpsFactorType::TYPE_A);
		SparseMatrixType Atranspose;
		transposeConjugate(Atranspose,A());
		assert(data_.size()==h.n_col());
		for (size_t b1=0;b1<data_.size();b1++)
			moveLeft(data_[b1],A,Atranspose,b1,h,dataPrev,symm);
	}

	void moveRight(const MpsFactorType& B,const MpoFactorType& h,const DataType& dataPrev,const SymmetryFactorType& symm)
	{
		assert(leftOrRight_ == PART_RIGHT);
		assert(B.type()==MpsFactorType::TYPE_B);
		SparseMatrixType Btranspose;
		transposeConjugate(Btranspose,B());
		if (h.n_row()!=data_.size())
			data_.resize(h.n_row());
		assert(h.n_row()==data_.size());
		for (size_t blm2=0;blm2<data_.size();blm2++)
			moveRight(data_[blm2],blm2,B,Btranspose,h,dataPrev,symm);
	}

	void moveRight(SparseMatrixType& m,
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
			moveRight(values,cols,alm2,blm2,B,Btranspose,h,dataPrev,symm);
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
		m.checkValidity();
	}

	void moveRight(std::vector<ComplexOrRealType>& values,
					 std::vector<size_t>& cols,
					 size_t alm2,
					 size_t blm2,
					 const MpsFactorType& B,
					 const SparseMatrixType& Btranspose,
					 const MpoFactorType& h,
					 const DataType& dataPrev,
					 const SymmetryFactorType& symm)
	{
		const SparseMatrixType& Bmatrix = B();
		assert(symm.right().split()==0 || symm.right().size()/symm.right().split()==dataPrev[0].row());
		assert(Btranspose.row()==symm.right().size());
		assert(dataPrev.size()<=h.n_col());
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
	size_t leftOrRight_;

}; // ContractedFactor

template<typename MatrixProductOperatorType>
std::ostream& operator<<(std::ostream& os,const ContractedFactor<MatrixProductOperatorType>& contractedFactor)
{
	os<<"SparseMatrices= "<<contractedFactor.data_.size()<<"\n";
	for (size_t i=0;i<contractedFactor.data_.size();i++) {
		typename MatrixProductOperatorType::MpsLocalType::MatrixType m;
		crsMatrixToFullMatrix(m,contractedFactor.data_[i]);
		std::cout<<m;
//		os<<contractedFactor.data_[i].row()<<"x"<<contractedFactor.data_[i].col()<<"    ";
	}
	os<<"\n";
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_FACTOR_H

