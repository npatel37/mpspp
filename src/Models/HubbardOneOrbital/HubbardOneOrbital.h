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
		 typename GeometryType,
		 typename ConcurrencyType>
class HubbardOneOrbital : public ModelBase<ParametersSolverType,
										   InputValidatorType,
										   GeometryType,
										   ConcurrencyType> {

	typedef ModelBase<ParametersSolverType,
					  InputValidatorType,
					  GeometryType,
					  ConcurrencyType> ModelBaseType;

	typedef typename ModelBaseType::MatrixProductOperatorType MatrixProductOperatorType;
	typedef typename ModelBaseType::SparseMatrixType SparseMatrixType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::VectorType VectorType;

public:

	HubbardOneOrbital(const ParametersSolverType& solverParams,
					  InputValidatorType& io,
					  const GeometryType& geometry,
					  ConcurrencyType& concurrency)
	: solverParams_(solverParams),
	  io_(io),
	  geometry_(geometry),
	  concurrency_(concurrency)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to set Hamiltonian here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	virtual const MatrixProductOperatorType& hamiltonian(size_t site) const
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

private:

	const ParametersSolverType& solverParams_;
	InputValidatorType& io_;
	const GeometryType& geometry_;
	ConcurrencyType& concurrency_;
	MatrixProductOperatorType hamiltonian_;

}; // HubbardOneOrbital

} // namespace Mpspp

/*@}*/
#endif // MODEL_HUBBARD_ONE_ORBITAL_H

