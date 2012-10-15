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

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT, PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:


	typedef typename ProgramGlobals::Vector<SparseMatrixType>::Type DataType;

	void init(const MpsFactorType& AorB,
			  const MpoFactorType& h,
			  size_t site,
			  size_t leftOrRight,
			  const ThisType* previous)
	{
		DataType prevfFirst;
		const DataType* prevf = 0;

		const SymmetryFactorType& symm = AorB.symm();

		if (previous==0) {
			SparseMatrixType m(1,1);
			m.makeDiagonal(1,1);
			prevfFirst.push_back(m);
			prevf = &prevfFirst;
		} else {
			prevf = &(previous->data_);
		}

		SparseMatrixType AorBtransp;
		transposeConjugate(AorBtransp,AorB());

		if (leftOrRight==PART_LEFT) {
			nextF(symm,AorB(),AorBtransp,h,*prevf,leftOrRight);
		} else {
			nextF(symm,AorBtransp,AorB(),h,*prevf,leftOrRight);
		}
	}

	//! From As (or Bs) and Ws reconstruct *this
	void update(const MpsFactorType& AorB,size_t direction)
	{
		if (direction==TO_THE_RIGHT) {
			updateLeft(AorB);
		} else {
			updateRight(AorB);
		}
	}

	const SparseMatrixType& operator()(size_t b1) const
	{
		return data_[b1];
	}

private:

	//! Eq.(197), page 63
	void nextF(const SymmetryFactorType& symm,const SparseMatrixType& A,const SparseMatrixType& Atransp,const MpoFactorType& h,const DataType& prevf,size_t leftOrRight)
	{
		size_t hilbertSize = h.hilbertSize();
		size_t leftBlockSize = A.row()/hilbertSize;

		data_.resize(prevf.size());
		DenseMatrixType matrix2(leftBlockSize,A.col());

		for (size_t bi=0;bi<prevf.size();bi++) {

			data_[bi].resize(A.col(),A.col());
			size_t counter = 0;
			size_t total = data_[bi].row();
			typename ProgramGlobals::Vector<int>::Type ptr(total,-1);
			typename ProgramGlobals::Vector<size_t>::Type index(total,0);
			VectorType temp(total,0.0);

			for (size_t ai=0;ai<A.col();ai++) {
				size_t itemp = 0;
				data_[bi].setRow(ai,counter);

				for (int k=Atransp.getRowPtr(ai);k<Atransp.getRowPtr(ai+1);k++) {
					size_t aim1si = Atransp.getCol(k);
					PairType aim1siP = symm.left().unpack(aim1si);
					size_t si = (leftBlockSize>1) ? aim1siP.second : aim1siP.first;

					middle197(matrix2,symm,si,bi,A,h,prevf,leftOrRight);
					size_t aim1 = (leftBlockSize>1) ? aim1siP.first : 0;
					for (size_t aip=0;aip<A.col();aip++) {
						ComplexOrRealType tmp = std::conj(Atransp.getValue(k)) * matrix2(aim1,aip);
						assert(aip<ptr.size());
						if (ptr[aip]<0) {
							ptr[aip]=itemp;
							temp[ptr[aip]]= tmp;
							index[ptr[aip]] = aip;
							itemp++;
						} else {
							temp[ptr[aip]]+=tmp;
						}
					}
				}

				for (size_t s=0;s<itemp;s++) {
					data_[bi].pushValue(temp[s]);
					data_[bi].pushCol(index[s]);
					ptr[index[s]] = -1;
				}
				counter += itemp;
				std::cerr<<"Testing ai="<<ai<<" out of "<<A.col()<<"\n";
			}

			data_[bi].setRow(total,counter);
			data_[bi].checkValidity();
		}
	}

	//! Eq.(197), page 63, middle bracket
	void middle197(DenseMatrixType& matrix2,const SymmetryFactorType& symm,size_t si,size_t bi,const SparseMatrixType& A,const MpoFactorType& h,const DataType& prevf,size_t leftOrRight) const
	{

		for (size_t i=0;i<matrix2.n_row();i++)
			for (size_t j=0;j<matrix2.n_col();j++)
				matrix2(i,j) = 0;

		for (size_t bim1=0;bim1<prevf.size();bim1++) {
			size_t bFirst = (leftOrRight==PART_LEFT) ? bim1  : bi;
			size_t bSecond =(leftOrRight==PART_LEFT) ? bi : bim1;
			const SparseMatrixType& wtmp = h(bFirst,bSecond);
			for (int ks=wtmp.getRowPtr(si);ks<wtmp.getRowPtr(si+1);ks++) {
				size_t sip = wtmp.getCol(ks);
				SparseMatrixType matrix(prevf[bim1].row(),A.col());
				inner197(matrix,symm,sip,A,prevf[bim1]);
				for (size_t aim1=0;aim1<matrix.row();aim1++) {
					for (int k=matrix.getRowPtr(aim1);k<matrix.getRowPtr(aim1+1);k++) {
						matrix2(aim1,matrix.getCol(k)) += wtmp.getValue(ks) * matrix.getValue(k);
					}
				}
			}
			//std::cerr<<"Testing bim1="<<bim1<<" out of "<<prevf.size()<<"\n";
		}
	}

	//! Eq.(197), page 63, inner bracket
	void inner197(SparseMatrixType& matrix,const SymmetryFactorType& symm,size_t sip,const SparseMatrixType& A,const SparseMatrixType& prevf) const
	{
		size_t counter = 0;
		//size_t total = matrix.row();
		typename ProgramGlobals::Vector<int>::Type ptr(A.col(),-1);
		typename ProgramGlobals::Vector<size_t>::Type index(matrix.col(),0);
		VectorType temp(matrix.col(),0.0);

		assert(prevf.row()==matrix.row());

		for (size_t aim1=0;aim1<prevf.row();aim1++) {
			size_t itemp = 0;
			matrix.setRow(aim1,counter);
			for (int k=prevf.getRowPtr(aim1);k<prevf.getRowPtr(aim1+1);k++) {
				size_t aipm1 = prevf.getCol(k);
				size_t aipm1sip = symm.left().pack(aipm1,sip);
				for (int k2=A.getRowPtr(aipm1sip);k2<A.getRowPtr(aipm1sip+1);k2++) {
					size_t aip = A.getCol(k2);
					ComplexOrRealType tmp = prevf.getValue(k) * A.getValue(k2);
					assert(aip<ptr.size());
					if (ptr[aip]<0) {
						ptr[aip]=itemp;
						temp[ptr[aip]]= tmp;
						index[ptr[aip]] = aip;
						itemp++;
					} else {
						temp[ptr[aip]]+=tmp;
					}
				}
			}
			for (size_t s=0;s<itemp;s++) {
				matrix.pushValue(temp[s]);
				matrix.pushCol(index[s]);
				ptr[index[s]] = -1;
			}
			counter += itemp;
		}
		matrix.setRow(matrix.row(),counter);
		matrix.checkValidity();
	}

	void updateLeft(const MpsFactorType& A)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to update(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	void updateRight(const MpsFactorType& B)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to update(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	DataType data_;

}; // ContractedFactor

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_FACTOR_H

