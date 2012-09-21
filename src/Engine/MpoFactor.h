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

#ifndef MPO_FACTOR_H
#define MPO_FACTOR_H

#include "ProgramGlobals.h"

namespace Mpspp {

template<typename RealType,typename ComplexOrRealType>
class MpoFactor {

	typedef typename ProgramGlobals::CrsMatrix<ComplexOrRealType>::Type SparseMatrixType;
	typedef typename ProgramGlobals::Matrix<SparseMatrixType>::Type MatrixType;

public:

	MpoFactor(size_t wdim1,size_t wdim2) : data_(wdim1,wdim2) {}

	MpoFactor(size_t wdim1) : data_(wdim1,1) {}

	const SparseMatrixType& operator()(size_t i,size_t j) const
	{
		assert(i<n_row() && j<n_col());
		return data_(i,j);
	}

	SparseMatrixType& operator()(size_t i,size_t j)
	{
		assert(i<n_row() && j<n_col());
		return data_(i,j);
	}

	const SparseMatrixType& operator()(size_t i) const
	{
		assert(i<n_row() && 0<n_col());
		return data_(i,0);
	}

	SparseMatrixType& operator()(size_t i)
	{
		assert(i<n_row() && 0<n_col());
		return data_(i,0);
	}

	size_t n_row() const { return data_.n_row(); }

	size_t n_col() const { return data_.n_col(); }

	size_t hilbertSize() const
	{
		assert(data_.n_row()>0 && data_.n_col()>0);
		return data_(0,0).row();
	}

private:

	MatrixType data_;

}; // MpoFactor

} // namespace Mpspp

/*@}*/
#endif // MPO_FACTOR_H
