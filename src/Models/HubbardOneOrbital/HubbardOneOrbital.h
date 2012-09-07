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

#ifndef MODEL_HUBBARD_ONE_ORBITAL_H
#define MODEL_HUBBARD_ONE_ORBITAL_H

#include "ModelBase.h"
#include "ParametersModelHubbard.h"
#include "ProgramGlobals.h"

/** Hamiltonian of the Hubbard Model:

  For a chain:
  left corner: [U_0 nup ndown, t_{01} c^\dagger,I]
  right corner: [I,c,U_{n-1} nup ndown]
  on middle site s: [ I                0          0 ]
					[ c                0          0 ]
					[ U_{s} nup ndown t_{s,s+1} c, I]

  For a ladder: ?

*/

namespace Mpspp {

template<typename ParametersSolverType,
		 typename InputValidatorType,
		 typename SymmetryLocalType,
		 typename GeometryType,
		 typename ConcurrencyType>
class HubbardOneOrbital : public ModelBase<ParametersSolverType,
										   InputValidatorType,
										   SymmetryLocalType,
										   GeometryType,
										   ConcurrencyType> {

	typedef ModelBase<ParametersSolverType,
					  InputValidatorType,
					  SymmetryLocalType,
					  GeometryType,
					  ConcurrencyType> ModelBaseType;

	typedef typename ModelBaseType::MatrixProductOperatorType MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::VectorType VectorType;
	typedef typename ModelBaseType::ComplexOrRealType ComplexOrRealType;
	typedef typename ParametersSolverType::RealType RealType;
	typedef typename ProgramGlobals::Matrix<ComplexOrRealType>::Type MatrixType;

	typedef ParametersModelHubbard<RealType> ParametersModelType;

public:

	enum {SPIN_UP,SPIN_DOWN};

	HubbardOneOrbital(const ParametersSolverType& solverParams,
					  InputValidatorType& io,
					  const GeometryType& geometry,
					  ConcurrencyType& concurrency)
	: solverParams_(solverParams),
	  io_(io),
	  geometry_(geometry),
	  concurrency_(concurrency),
	  mp_(io),
	  hamiltonian_(geometry_.numberOfSites())
	{
		// FIXME: CONNECT WITH THE GEOMETRY HERE!!
		RealType tiip1 = 1.0;
		size_t n = hamiltonian_.size();
		size_t wdim = 4;
		size_t hilbert = 4;

		SparseMatrixType identity(hilbert,hilbert);
		identity.makeDiagonal(hilbert,1.0);
		SparseMatrixType cup(hilbert,hilbert);
		fillDestructionMatrix(cup,SPIN_UP);
		SparseMatrixType cdaggerUp(hilbert,hilbert);
		transposeConjugate(cdaggerUp,cup);

		SparseMatrixType cdown(hilbert,hilbert);
		fillDestructionMatrix(cdown,SPIN_DOWN);
		SparseMatrixType cdaggerDown(hilbert,hilbert);
		transposeConjugate(cdaggerDown,cdown);

		SparseMatrixType nupndown = (cdaggerUp*cup) * (cdaggerDown*cdown);

		MpoFactorType mleft(wdim);
		mleft(0) = mp_.hubbardU[0]* nupndown;
		mleft(1) = tiip1 * cdaggerUp;
		mleft(2) = tiip1 * cdaggerDown;
		mleft(3) = identity;
		hamiltonian_(0)=mleft;

		for (size_t i=1;i<n-1;i++) {
			MpoFactorType m(wdim,wdim);
			m(0,0) = identity;
			m(1,0) = cup;
			m(2,0) = cdown;
			m(3,0) = mp_.hubbardU[i]* nupndown;
			m(3,1) = tiip1 * cdaggerUp;
			m(3,2) = tiip1 * cdaggerDown;
			m(3,2) = identity;
			hamiltonian_(i)=m;
		}

		MpoFactorType mright(wdim);
		mright(3) = mp_.hubbardU[n-1]* nupndown;
		mright(1) = cup;
		mright(2) = cdown;
		mright(0) = identity;
		hamiltonian_(n-1)=mright;
	}

	virtual const MatrixProductOperatorType& hamiltonian() const
	{
		return hamiltonian_;
	}

	virtual void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need fullHamiltonian here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	virtual void matrixVectorProduct(VectorType& x,const VectorType& y,const ModelHelperType& modelHelper) const
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need matrixVectorProduct here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	virtual const ParametersSolverType& solverParams() const { return solverParams_; }

	virtual const GeometryType& geometry() const { return geometry_; }

private:

	void fillDestructionMatrix(SparseMatrixType& cm,size_t spin) const
	{
		MatrixType m(4,4);
		if (spin==SPIN_UP) {
			m(1,0) = 1;
			m(3,2) = 1;
		} else {
			m(2,0) = 1;
			m(3,1) = -1;
		}
		fullMatrixToCrsMatrix(cm,m);
	}

	const ParametersSolverType& solverParams_;
	InputValidatorType& io_;
	const GeometryType& geometry_;
	ConcurrencyType& concurrency_;
	ParametersModelType mp_;
	MatrixProductOperatorType hamiltonian_;

}; // HubbardOneOrbital

} // namespace Mpspp

/*@}*/
#endif // MODEL_HUBBARD_ONE_ORBITAL_H

