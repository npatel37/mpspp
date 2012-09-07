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
	typedef typename MatrixProductStateType::MpsFactorType  MpsFactorType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MatrixProductStateType::ComplexOrRealType ComplexOrRealType;
	typedef ContractedFactor<MatrixProductOperatorType> ThisType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

public:

	typedef typename ProgramGlobals::CrsMatrix<ComplexOrRealType>::Type SparseMatrixType;

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

	typename ProgramGlobals::Vector<SparseMatrixType>::Type data_;

}; // ContractedFactor

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_FACTOR_H

