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

public:


	typedef typename ProgramGlobals::Vector<SparseMatrixType>::Type DataType;

	void init(const MpsFactorType& AorB,
			  const MpoFactorType& h,
			  size_t site,
			  size_t leftOrRight,
			  const ThisType* previous)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to set data_ here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
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
	void nextF(DataType& thisf,const MpsFactorType& A,const MpsFactorType& Adagger,const MpoFactorType& h,const DataType& prevf) const
	{
		size_t hilbertSize = h.hilbertSize();
		size_t leftBlockSize = A.row()/hilbertSize;

		for (size_t bi=0;bi<prevf.size();bi++) {

			size_t counter = 0;
			size_t total = thisf[bi].row();
			typename ProgramGlobals::Vector<int>::Type ptr(total,-1);
			typename ProgramGlobals::Vector<size_t>::Type index(total,0);
			VectorType temp(total,0.0);
			for (size_t ai=0;ai<leftBlockSize;ai++) {
				size_t itemp = 0;
				thisf[bi].setRow(ai,counter);

				for (int k=Adagger.getRowPtr(ai);k<Adagger.getRowPtr(ai+1);k++) {
					size_t aim1si = Adagger.getCol(k);
					PairType aim1siP = Adagger.symm().left().unpack(aim1si);
					size_t si = aim1siP.second;
					DenseMatrixType matrix2(leftBlockSize,A.col());
					middle197(matrix2,si,bi,A,prevf);
					size_t aim1 = aim1siP.first;
					for (size_t aip=0;aip<A.getCol();aip++) {
						ComplexOrRealType tmp = Adagger.getValue(k) * matrix2(aim1,aip);
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
					thisf[bi].pushValue(temp[s]);
					thisf[bi].pushCol(index[s]);
					ptr[index[s]] = -1;
				}
				counter += itemp;
			}

			thisf[bi].setRow(total,counter);
			thisf[bi].checkValidity();
		}
	}

	//! Eq.(197), page 63, middle bracket
	void middle197(DenseMatrixType& matrix2,size_t si,size_t bi,const SparseMatrixType& A,const MpoFactorType& h,const DataType& prevf) const
	{
		for (size_t bim1=0;bim1<prevf.size();bim1++) {
			const SparseMatrixType& wtmp = h(bim1,bi);
			for (size_t ks=wtmp.getRowPtr(si);ks<wtmp.getRowPtr(si+1);ks++) {
				size_t sip = wtmp.getCol(ks);
				SparseMatrixType matrix(prevf.row(),A.col());
				inner197(matrix,sip,A,prevf[bim1]);
				for (size_t aim1=0;aim1<matrix.row();aim1++) {
					for (size_t k=matrix.getRowPtr(aim1);k<matrix.getRowPtr(aim1+1);k++) {
						matrix2(aim1,matrix.getCol(k)) += wtmp.getValue(ks) * matrix.getValue(k);
					}

				}
			}
		}
	}

	//! Eq.(197), page 63, inner bracket
	void inner197(SparseMatrixType& matrix,size_t sip,const SparseMatrixType& A,const SparseMatrixType& prevf) const
	{
		size_t counter = 0;
		size_t total = matrix.row();
		typename ProgramGlobals::Vector<int>::Type ptr(total,-1);
		typename ProgramGlobals::Vector<size_t>::Type index(total,0);
		VectorType temp(total,0.0);
		for (size_t aim1=0;aim1<prevf.row();aim1++) {
			size_t itemp = 0;
			matrix.setRow(aim1,counter);
			for (size_t k=prevf.getRowPtr(aim1);k<prevf.getRowPtr(aim1+1);k++) {
				size_t aipm1 = prevf.getCol(k);
				size_t aipm1sip = A.symm().left().pack(aipm1,sip);
				for (size_t k2=A.getRow(aipm1sip);k2<A.getRow(aipm1sip+1);k2++) {
					size_t aip = A.getCol(k2);
					ComplexOrRealType tmp = prevf.getValue(k) * A.getValue(k2);
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
		matrix.setRow(total,counter);
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

