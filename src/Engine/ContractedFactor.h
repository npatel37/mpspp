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

//	ContractedFactor()
//	{
//		SparseMatrixType m(1,1);
//		m.makeDiagonal(1,1);
//		data_.push_back(m);
//	}

	ContractedFactor(const MpsFactorType& AorB,const MpoFactorType& h,size_t site,size_t leftOrRight,ThisType* dataPrev)
		: data_(h.n_col()),site_(site),leftOrRight_(leftOrRight)
	{
		SparseMatrixType Atranspose;
		transposeConjugate(Atranspose,AorB());
		for (size_t b1=0;b1<data_.size();b1++) {
			if (leftOrRight == PART_RIGHT) {
				initRight2(data_[b1],AorB,b1,h);
			} else {
				initLeft2(data_[b1],AorB,Atranspose,b1,h,dataPrev);
			}
		}
	}

//	void initLeft(const MpsFactorType& AorB,
//			  const MpoFactorType& h,
//			  size_t site,
//			  const ThisType* previous)
//	{
//		//DataType prevfFirst;
//		const DataType* prevf = 0;

//		const SymmetryFactorType& symm = AorB.symm();

//		if (previous==0) {
//			SparseMatrixType m(1,1);
//			m.makeDiagonal(1,1);
//			data_.push_back(m);
//			return;
//			//prevfFirst.push_back(m);
//			//prevf = &prevfFirst;
//		} else {
//			prevf = &(previous->data_);
//		}

//		SparseMatrixType AorBtransp;
//		transposeConjugate(AorBtransp,AorB());

//		//if (leftOrRight==PART_LEFT) {
//			nextF(symm,AorB(),AorBtransp,h,*prevf,PART_LEFT);
////		} else {
////			nextF(symm,AorBtransp,AorB(),h,*prevf,leftOrRight);
////		}
//	}

	//! From As (or Bs) and Ws reconstruct *this
	void update(const MpsFactorType& AorB)
	{
		if (AorB.type()==MpsFactorType::TYPE_A) {
			updateLeft(AorB);
		} else {
			updateRight(AorB);
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

	void initLeft2(SparseMatrixType& m,const MpsFactorType& AorB,const SparseMatrixType& Atranspose,size_t b,const MpoFactorType& h,ThisType* dataPrev)
	{
		if (site_==0) {
			contractedFactor0(m,AorB,b,h);
			return;
		}

		std::cerr<<"Start initLeft2\n";

		assert(dataPrev!=0);

		const SymmetryFactorType& symm = AorB.symm();
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

//	//! Eq.(197), page 63
//	void nextF(const SymmetryFactorType& symm,const SparseMatrixType& A,const SparseMatrixType& Atransp,const MpoFactorType& h,const DataType& prevf,size_t leftOrRight)
//	{
//		size_t hilbertSize = h.hilbertSize();
//		size_t leftBlockSize = A.row()/hilbertSize;
//		size_t hsize = h(0,0).col();
//		size_t hsize2 = h.n_row();

//		data_.resize(hsize2);

//		DenseMatrixType dataMatrix(A.col(),A.col());

//		auxStorage_.resize(hsize,hsize2);
//		auxStorage2_.resize(hsize,hsize2);

//		size_t siSize = (leftBlockSize>1) ? hilbertSize : symm.left().split();

//		for (size_t bi=0;bi<hsize2;bi++)
//			for (size_t si=0;si<siSize;si++)
//				saveMiddle197(symm,si,bi,A,h,prevf,leftOrRight);


//		for (size_t bi=0;bi<hsize2;bi++) {
//			dataMatrix.setTo(0);
//			for (size_t ai=0;ai<A.col();ai++) {
//				for (int k=Atransp.getRowPtr(ai);k<Atransp.getRowPtr(ai+1);k++) {
//					size_t aim1si = Atransp.getCol(k);
//					PairType aim1siP = symm.left().unpack(aim1si);
//					size_t si = (leftBlockSize>1) ? aim1siP.second : aim1siP.first;

//					//middle197(matrix2,symm,si,bi,A,h,prevf,leftOrRight);
//					const DenseMatrixType& matrix2 = auxStorage2_(si,bi);
//					size_t aim1 = (leftBlockSize>1) ? aim1siP.first : 0;
//					for (size_t aip=0;aip<A.col();aip++) {
//						dataMatrix(ai,aip) += std::conj(Atransp.getValue(k)) * matrix2(aim1,aip);
//					}
//				}
//			}
//			fullMatrixToCrsMatrix(data_[bi],dataMatrix);
//			std::cerr<<"nonZero for bi="<<bi<<" is "<<data_[bi].nonZero()<<" row="<<data_[bi].row()<<"\n";
//		}

//		auxStorage_.reset(0,0);
//		auxStorage2_.reset(0,0);
//	}

//	void saveMiddle197(const SymmetryFactorType& symm,size_t si,size_t bi,const SparseMatrixType& A,const MpoFactorType& h,const DataType& prevf,size_t leftOrRight)
//	{
//		if (auxStorage2_(si,bi).n_row()!=0) return;

//		size_t hilbertSize = h.hilbertSize();
//		size_t leftBlockSize = A.row()/hilbertSize;

//		DenseMatrixType matrix2(leftBlockSize,A.col());

//		middle197(matrix2,symm,si,bi,A,h,prevf,leftOrRight);
////		SparseMatrixType matrixSparse(matrix2);
//		auxStorage2_(si,bi) = matrix2;
//	}

//	//! Eq.(197), page 63, middle bracket
//	void middle197(DenseMatrixType& matrix2,const SymmetryFactorType& symm,size_t si,size_t bi,const SparseMatrixType& A,const MpoFactorType& h,const DataType& prevf,size_t leftOrRight)
//	{
//		matrix2.setTo(0);

//		size_t hsize2 = (leftOrRight==PART_LEFT) ? h.n_row() : h.n_col();

//		for (size_t bim1=0;bim1<hsize2;bim1++) {
//			size_t bFirst = (leftOrRight==PART_LEFT) ? bim1  : bi;
//			size_t bSecond =(leftOrRight==PART_LEFT) ? bi : bim1;
//			const SparseMatrixType& wtmp = h(bFirst,bSecond);
//			if (wtmp.row()==0) continue;

//			for (int ks=wtmp.getRowPtr(si);ks<wtmp.getRowPtr(si+1);ks++) {
//				size_t sip = wtmp.getCol(ks);
//				saveInner197(symm,sip,bim1,A,prevf);
//				//SparseMatrixType matrixSparse(matrix);
//				const SparseMatrixType& matrixSparse = auxStorage_(sip,bim1);
//				for (size_t aim1=0;aim1<matrixSparse.row();aim1++) {
//					size_t start = matrixSparse.getRowPtr(aim1);
//					size_t end = matrixSparse.getRowPtr(aim1+1);
//					for (size_t k=start;k<end;++k) {
//						matrix2(aim1,matrixSparse.getCol(k)) += wtmp.getValue(ks) * matrixSparse.getValue(k);
//					}
//				}
//			}
//			//std::cerr<<"Testing bim1="<<bim1<<" out of "<<hsize2<<"\n";
//		}
//	}

//	//! Eq.(197), page 63, inner bracket
//	//! Takes into account truncation
//	void inner197(DenseMatrixType& matrixTemp,const SymmetryFactorType& symm,size_t sip,const SparseMatrixType& A,const SparseMatrixType& prevf) const
//	{
//		size_t leftSize = std::min(prevf.row(),symm.left().split());

//		for (size_t aim1=0;aim1<leftSize;aim1++) {
//			for (size_t jj=0;jj<matrixTemp.n_col();++jj) matrixTemp(aim1,jj)=0;
//			for (int k=prevf.getRowPtr(aim1);k<prevf.getRowPtr(aim1+1);k++) {
//				size_t aipm1 = prevf.getCol(k);
//				if (aipm1>=leftSize) continue;
//				size_t aipm1sip = symm.left().pack(aipm1,sip);
//				for (int k2=A.getRowPtr(aipm1sip);k2<A.getRowPtr(aipm1sip+1);k2++) {
//					size_t aip = A.getCol(k2);
//					matrixTemp(aim1,aip) += prevf.getValue(k) * A.getValue(k2);
//				}
//			}
//		}
//	}

//	void saveInner197(const SymmetryFactorType& symm,size_t sip,size_t bim1,const SparseMatrixType& A,const DataType& prevf)
//	{
//		if (auxStorage_(sip,bim1).row()!=0) return;

////		assert(hsize2>0);
////		size_t prevfRow = prevf[0].row();
////		for (size_t bim2=0;bim2<hsize2;bim2++) {
////			assert(prevf[bim2].row()==prevfRow);
////		}

//		DenseMatrixType matrixTemp(symm.left().split(),A.col());

//		inner197(matrixTemp,symm,sip,A,prevf[bim1]);
//		SparseMatrixType matrixSparse(matrixTemp);
//		auxStorage_(sip,bim1) = matrixSparse;
//	}

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
	size_t site_;
	size_t leftOrRight_;
//	typename ProgramGlobals::Matrix<SparseMatrixType>::Type auxStorage_;
//	typename ProgramGlobals::Matrix<DenseMatrixType>::Type auxStorage2_;

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

