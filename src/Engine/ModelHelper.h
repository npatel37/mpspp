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

#ifndef MODEL_HELPER_H
#define MODEL_HELPER_H
#include "ProgramGlobals.h"

namespace Mpspp {

template<typename LeftRightSuperType>
class ModelHelper {

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

public:

	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename LeftRightSuperType::ComplexOrRealType ComplexOrRealType;
	typedef typename LeftRightSuperType::ContractedPartType ContractedPartType;
	typedef typename ContractedPartType::ContractedFactorType ContractedFactorType;
	typedef typename LeftRightSuperType::MatrixProductStateType MatrixProductStateType;
	typedef typename MatrixProductStateType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef typename ContractedPartType::SparseMatrixType SparseMatrixType;
	typedef typename ProgramGlobals::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename LeftRightSuperType::MatrixProductOperatorType MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;

	ModelHelper(const LeftRightSuperType& lrs,
				size_t symmetrySector,
				size_t currentSite,
				size_t direction,
				const MpoFactorType& hamiltonian)
	: lrs_(lrs),
	  symmetrySector_(symmetrySector),
	  currentSite_(currentSite),
	  direction_(direction),
	  hamiltonian_(hamiltonian),
	  symmetry_(lrs_.symmetry()(currentSite_))
	{}

	size_t size() const
	{
		return lrs_.symmetry()(currentSite_).super().partitionSize(symmetrySector_);
	}

	size_t symmetrySector() const { return symmetrySector_; }

	size_t hilbertSize() const { return hamiltonian_(0,0).row(); }

	const LeftRightSuperType& lrs() const { return lrs_; }

	const MpoFactorType& hamiltonian() const { return hamiltonian_; }

	const SymmetryFactorType& symmetry() const { return symmetry_; }

//	const ContractedFactorType& contractedFactorLeft() const
//	{
//		return lrs_.contracted()(currentSite_,ProgramGlobals::PART_LEFT);
//	}

//	const ContractedFactorType& contractedFactorRight() const
//	{
//		return lrs_.contracted()(currentSite_,ProgramGlobals::PART_RIGHT);
//	}

	//! Eq. (201) but very modified
	void matrixVectorProduct(VectorType& x,const VectorType& y) const
	{
		size_t offset = symmetry_.super().partitionOffset(symmetrySector_);
		size_t total = symmetry_.super().partitionSize(symmetrySector_);
//		size_t nright = symmetry_.right().block().size();
//		assert(nright>0);
//		size_t center = symmetry_.right().block()[nright-1];

		if (currentSite_==0) return matrixVectorProduct0(x,y);

		const ContractedFactorType& cL = lrs_.contracted()(currentSite_-1,ProgramGlobals::PART_LEFT);
		const ContractedFactorType& cR = lrs_.contracted()(currentSite_,ProgramGlobals::PART_RIGHT);
		const SymmetryFactorType& symm = symmetry_;
		for (size_t blm1=0;blm1<cL.size();blm1++) {
			const SparseMatrixType& l1 = cL(blm1);
			for (size_t bl=0;bl<cR.size();bl++) {
				const SparseMatrixType& w =  hamiltonian_(blm1,bl);
				if (w.row()==0) continue;
//				SparseMatrixType w;
//				transposeConjugate(w,w1);
				const SparseMatrixType& r1 = cR(bl);
				for (size_t i=0;i<total;i++) {
					PairType ab = symm.super().unpack(i+offset);
					size_t alm1=0;
					size_t sigmaL=0;
					size_t alB=0;
					if (direction_==TO_THE_RIGHT) {
						PairType tmpPair1 = symm.left().unpack(ab.first);
						alm1=tmpPair1.first;
						sigmaL=tmpPair1.second;
						alB = ab.second;
					} else {
						PairType tmpPair1 = symm.right().unpack(ab.second);
						alB=tmpPair1.second;
						sigmaL=tmpPair1.first;
						alm1=ab.first;
					}

					for (int k1=l1.getRowPtr(alm1);k1<l1.getRowPtr(alm1+1);k1++) {
						size_t alm1p=l1.getCol(k1);
						for (int kw=w.getRowPtr(sigmaL);kw<w.getRowPtr(sigmaL+1);kw++) {
							size_t sigmaLp=w.getCol(kw);
							for (int k2=r1.getRowPtr(alB);k2<r1.getRowPtr(alB+1);k2++) {
								size_t alBp = r1.getCol(k2);
								size_t j = 0;
								if (direction_==TO_THE_RIGHT) {
									size_t tmp1 = symm.left().pack(alm1p,sigmaLp);
									j = symm.super().pack(tmp1,alBp);
								} else {
									size_t tmp1 = symm.right().pack(sigmaLp,alBp);
									j = symm.super().pack(alm1,tmp1);
								}
								if (j<offset || j>=offset+total) continue;
								x[i] += y[j-offset]*l1.getValue(k1)*w.getValue(kw)*r1.getValue(k2);
							} // k2 right
						} // kw Hamiltonian
					} // k1 left
				} // symmetry sector
			} // bl
		} // blm1
	}

	//! Used only for stored option
	void fullHamiltonian(SparseMatrixType& matrix) const
	{
		size_t offset = symmetry_.super().partitionOffset(symmetrySector_);
		size_t total = symmetry_.super().partitionSize(symmetrySector_);
//		size_t nright = symmetry_.right().block().size();
//		assert(nright>0);
//		size_t center = symmetry_.right().block()[nright-1];

		if (currentSite_==0) return fullHamiltonian0(matrix);

		matrix.resize(total,total);
		const ContractedFactorType& cL = lrs_.contracted()(currentSite_-1,ProgramGlobals::PART_LEFT);
		const ContractedFactorType& cR = lrs_.contracted()(currentSite_,ProgramGlobals::PART_RIGHT);
		const SymmetryFactorType& symm = symmetry_;
		VectorType v(total,0);
		size_t counter = 0;
		for (size_t i=0;i<total;i++) {
			matrix.setRow(i,counter);
			for (size_t blm1=0;blm1<cL.size();blm1++) {
				const SparseMatrixType& l1 = cL(blm1);
				for (size_t bl=0;bl<cR.size();bl++) {
					const SparseMatrixType& w = hamiltonian_(blm1,bl);
					if (w.row()==0) continue;
//					SparseMatrixType w;
//					transposeConjugate(w,w1);
					const SparseMatrixType& r1 = cR(bl);

					PairType ab = symm.super().unpack(i+offset);
					size_t alm1=0;
					size_t sigmaL=0;
					size_t alB=0;
					if (direction_==TO_THE_RIGHT) {
						PairType tmpPair1 = symm.left().unpack(ab.first);
						alm1=tmpPair1.first;
						sigmaL=tmpPair1.second;
						alB = ab.second;
					} else {
						PairType tmpPair1 = symm.right().unpack(ab.second);
						alB=tmpPair1.second;
						sigmaL=tmpPair1.first;
						alm1=ab.first;
					}

					for (int k1=l1.getRowPtr(alm1);k1<l1.getRowPtr(alm1+1);k1++) {
						size_t alm1p=l1.getCol(k1);
						for (int kw=w.getRowPtr(sigmaL);kw<w.getRowPtr(sigmaL+1);kw++) {
							size_t sigmaLp=w.getCol(kw);
							for (int k2=r1.getRowPtr(alB);k2<r1.getRowPtr(alB+1);k2++) {
								size_t alBp = r1.getCol(k2);
								size_t j = 0;
								if (direction_==TO_THE_RIGHT) {
									size_t tmp1 =  symm.left().pack(alm1p,sigmaLp);
									j = symm.super().pack(tmp1,alBp);
								} else {
									size_t tmp1 = symm.right().pack(sigmaLp,alBp);
									j = symm.super().pack(alm1,tmp1);
								}
								if (j<offset || j>=offset+total) continue;
								v[j-offset] += l1.getValue(k1)*w.getValue(kw)*r1.getValue(k2);
							} // k2 right
						} // kw Hamiltonian
					} // k1 left
				} // bl
			} // blm1
			for (size_t j=0;j<v.size();j++) {
				if (fabs(v[j])<1e-6) continue;
				matrix.pushCol(j);
				matrix.pushValue(v[j]);
				v[j]=0;
				counter++;
			}
		} // symmetry sector
		matrix.setRow(total,counter);
		matrix.checkValidity();
	}

private:

	//! Eq. (201) but very modified
	void matrixVectorProduct0(VectorType& x,const VectorType& y) const
	{
		size_t offset = symmetry_.super().partitionOffset(symmetrySector_);
		size_t total = symmetry_.super().partitionSize(symmetrySector_);

		const ContractedFactorType& cR = lrs_.contracted()(currentSite_,ProgramGlobals::PART_RIGHT);
		const SymmetryFactorType& symm = symmetry_;

		for (size_t bl=0;bl<cR.size();bl++) {
			const SparseMatrixType& w =  hamiltonian_(bl);
			if (w.row()==0) continue;
//			SparseMatrixType w;
//			transposeConjugate(w,w1);
			const SparseMatrixType& r1 = cR(bl);
			for (size_t i=0;i<total;i++) {
				PairType ab = symm.super().unpack(i+offset);
				size_t a1 = ab.first;
				size_t sigma1 = ab.second;

				for (int kw=w.getRowPtr(sigma1);kw<w.getRowPtr(sigma1+1);kw++) {
					size_t sigma1p=w.getCol(kw);
					for (int k2=r1.getRowPtr(a1);k2<r1.getRowPtr(a1+1);k2++) {
						size_t a1p = r1.getCol(k2);
						size_t j = symm.super().pack(a1p,sigma1p);
						if (j<offset || j>=offset+total) continue;
						x[i] += y[j-offset]*w.getValue(kw)*r1.getValue(k2);
					} // k2 right
				} // kw Hamiltonian
			} // symmetry sector
		} // bl
	}

	//! Used only for stored option
	void fullHamiltonian0(SparseMatrixType& matrix) const
	{
		size_t offset = symmetry_.super().partitionOffset(symmetrySector_);
		size_t total = symmetry_.super().partitionSize(symmetrySector_);

		matrix.resize(total,total);
		const ContractedFactorType& cR = lrs_.contracted()(currentSite_,ProgramGlobals::PART_RIGHT);
		const SymmetryFactorType& symm = symmetry_;
		VectorType v(total,0);
		size_t counter = 0;
		for (size_t i=0;i<total;i++) {
			matrix.setRow(i,counter);
			for (size_t bl=0;bl<cR.size();bl++) {
				const SparseMatrixType& w = hamiltonian_(bl);
				if (w.row()==0) continue;
//				SparseMatrixType w;
//				transposeConjugate(w,w1);
				const SparseMatrixType& r1 = cR(bl);

				PairType ab = symm.super().unpack(i+offset);
				size_t a1 = ab.first;
				size_t sigma1 = ab.second;

				for (int kw=w.getRowPtr(sigma1);kw<w.getRowPtr(sigma1+1);kw++) {
					size_t sigma1p=w.getCol(kw);
					for (int k2=r1.getRowPtr(a1);k2<r1.getRowPtr(a1+1);k2++) {
						size_t a1p = r1.getCol(k2);
						size_t j = symm.super().pack(a1p,sigma1p);
						if (j<offset || j>=offset+total) continue;
						v[j-offset] += w.getValue(kw)*r1.getValue(k2);
					} // k2 right
				} // kw Hamiltonian
			} // bl

			for (size_t j=0;j<v.size();j++) {
				if (fabs(v[j])<1e-6) continue;
				matrix.pushCol(j);
				matrix.pushValue(v[j]);
				v[j]=0;
				counter++;
			}
		} // symmetry sector
		matrix.setRow(total,counter);
		matrix.checkValidity();
	}


	const LeftRightSuperType& lrs_;
	size_t symmetrySector_;
	size_t currentSite_;
	size_t direction_;
	size_t hilbertSize_;
	const MpoFactorType& hamiltonian_;
	const SymmetryFactorType& symmetry_;

}; // ModelHelper

} // namespace Mpspp

/*@}*/
#endif // MODEL_HELPER_H

