/*
Copyright (c) 2012-2013, UT-Battelle, LLC
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

#ifndef SYMMETRY_HELPER_H
#define SYMMETRY_HELPER_H

#include "ProgramGlobals.h"
#include "MpsLocal.h"
#include "MpoFactor.h"

namespace Mpspp {

template<typename SymmetryFactorType>
class SymmetryHelper {

	typedef ProgramGlobals::Vector<size_t>::Type VectorIntegerType;

	static VectorIntegerType electronsFromQn_;
	static const int MAX_SITES = ProgramGlobals::MAX_SITES;

public:

	template<typename SomeFermionSignType>
	SymmetryHelper(const SomeFermionSignType& fermionSign,const SymmetryFactorType& symm)
		: symm_(symm)
	{
		VectorIntegerType qn = fermionSign.quantumNumbers();
		electronsOneSite_.resize(qn.size());
		for (size_t i=0;i<qn.size();i++)
			electronsOneSite_[i] = fermionSign.electronsFromQn(qn[i]);

		if (electronsFromQn_.size()>0) return;

		electronsFromQn_.resize(MAX_SITES*MAX_SITES);

		for (size_t q=0;q<electronsFromQn_.size();q++)
			electronsFromQn_[q] = fermionSign.electronsFromQn(q);
	}

	const SymmetryFactorType& symm() const
	{
		return symm_;
	}

	size_t electronsFromState(size_t state) const
	{
		return electronsOneSite_[state];
	}

	size_t electronsFromQn(size_t q) const
	{
		return electronsFromQn_[q];
	}

private:

	const SymmetryFactorType& symm_;
	VectorIntegerType electronsOneSite_;

}; // SymmetryHelper

template<typename SymmetryFactorType>
typename SymmetryHelper<SymmetryFactorType>::VectorIntegerType SymmetryHelper<SymmetryFactorType>::electronsFromQn_;

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_HELPER_H

