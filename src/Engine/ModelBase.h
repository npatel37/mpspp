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

#ifndef MODEL_BASE_H
#define MODEL_BASE_H
#include "MatrixProductOperator.h"
#include "ModelHelper.h"
#include "ReflectionSymmetryEmpty.h"
#include "ContractedPart.h"

namespace Mpspp {

template<typename ParametersSolverType_,
		 typename InputValidatorType_,
		 typename SymmetryLocalType_,
		 typename GeometryType_,
		 typename ConcurrencyType_>
class ModelBase {

public:

	typedef ParametersSolverType_ ParametersSolverType;
	typedef typename ParametersSolverType::ComplexOrRealType ComplexOrRealType;
	typedef InputValidatorType_ InputValidatorType;
	typedef SymmetryLocalType_ SymmetryLocalType;
	typedef GeometryType_ GeometryType;
	typedef ConcurrencyType_ ConcurrencyType;
	typedef MatrixProductOperator<ComplexOrRealType,SymmetryLocalType> MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ProgramGlobals::Vector<RealType>::Type VectorType;
	typedef ContractedPart<MatrixProductOperatorType> ContractedPartType;
	typedef ModelHelper<ContractedPartType> ModelHelperType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef ReflectionSymmetryEmpty<SparseMatrixType> ReflectionSymmetryType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type MatrixType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;

	virtual ~ModelBase() {}

	virtual const MatrixProductOperatorType& hamiltonian() const=0;

	virtual const ParametersSolverType& solverParams() const=0;

	virtual const GeometryType& geometry() const=0;

//	virtual void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
//	{
//		typedef typename SymmetryLocalType::PairType PairType;
//		size_t symmetrySector = modelHelper.symmetrySector();

//		const SymmetryFactorType& symm = modelHelper.symmetry();

//		size_t total = symm.super().partitionSize(symmetrySector);
//		size_t offset = symm.super().partitionOffset(symmetrySector);

//		const MpoFactorType& hamiltonian = modelHelper.hamiltonian();

//		matrix.resize(total,total);
//		size_t counter = 0;
//		typename ProgramGlobals::Vector<int>::Type ptr(total,-1);
//		typename ProgramGlobals::Vector<size_t>::Type index(total,0);
//		VectorType temp(total,0.0);

//		for (size_t i=0;i<total;i++) {
//			matrix.setRow(i,counter);
//			size_t itemp = 0;
//			PairType iLa2 = symm.super().unpack(i+offset);
//			size_t iL = iLa2.first;
//			size_t a2 = iLa2.second;
//			PairType a1sigma2 = symm.left().unpack(iL);
//			size_t sigma2 = a1sigma2.second;
//			size_t a1 = a1sigma2.first;
//			for (size_t b1=0;b1<hamiltonian.n_row();b1++) {
//				const SparseMatrixType& cLm = modelHelper.contractedFactorLeft()(b1);
//				for (int k1=cLm.getRowPtr(a1);k1<cLm.getRowPtr(a1+1);k1++) {
//					size_t a1prime = cLm.getCol(k1);
//					for (size_t b2=0;b2<hamiltonian.n_col();b2++) {
//						const MatrixType& wm = hamiltonian(b1,b2);
//						const SparseMatrixType& cRm = modelHelper.contractedFactorRight()(b2);
//						for (int k2=cRm.getRowPtr(a2);k2<cRm.getRowPtr(a2+1);k2++) {
//							size_t a2prime = cLm.getCol(k2);
//							for (size_t sigma2prime = 0;sigma2prime<modelHelper.hilbertSize();sigma2prime++) {
//								size_t iLprime = symm.left().pack(a1prime,sigma2prime);
//								size_t iprime = symm.super().pack(iLprime,a2prime);
//								assert(iprime>=offset && iprime<total+offset);
//								iprime -= offset;
//								assert(iprime>=0);
//								ComplexOrRealType tmp = cLm.getValue(k1) * wm(sigma2,sigma2prime) * cRm.getValue(k2);
//								if (ptr[iprime]<0) {
//									ptr[iprime]=itemp;
//									temp[ptr[iprime]]= tmp;
//									index[ptr[iprime]] = iprime;
//									itemp++;
//								} else {
//									temp[ptr[iprime]]+=tmp;
//								}
//							}
//						}
//					}
//				}
//			}
//			for (size_t s=0;s<itemp;s++) {
//				matrix.pushValue(temp[s]);
//				matrix.pushCol(index[s]);
//				ptr[index[s]] = -1;
//			}
//			counter += itemp;
//		}
//		matrix.setRow(total,counter);
//		matrix.checkValidity();
//	}


	virtual void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
	{
		modelHelper.fullHamiltonian(matrix);
	}

	virtual void matrixVectorProduct(VectorType& x,const VectorType& y,const ModelHelperType& modelHelper) const
	{
		modelHelper.matrixVectorProduct(x,y);
	}

	virtual void getOneSite(std::vector<size_t>& quantumNumbers,size_t site) const=0;

}; // ModelBase

} // namespace Mpspp

/*@}*/
#endif // MODEL_BASE_H

