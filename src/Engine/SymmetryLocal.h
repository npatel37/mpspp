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

#ifndef SYMMETRY_LOCAL_H
#define SYMMETRY_LOCAL_H

#include "ProgramGlobals.h"
#include "SymmetryFactor.h"

namespace Mpspp {

class SymmetryLocal {

public:

	typedef SymmetryFactor SymmetryFactorType;
	typedef typename SymmetryFactorType::PairType PairType;
	typedef typename SymmetryFactorType::IoInputType IoInputType;

	SymmetryLocal(IoInputType& io)
	{
		size_t n = 0;
		io.readline(n,"TotalNumberOfSites=");
		size_t nk = 0;
		io.readline(nk,"HilbertOneSite=");
		assert(n>2);
		for (size_t i=0;i<n-1;i++) {
			SymmetryFactorType f(io,nk);
			data_.push_back(f);
//			if (i==0) data_.push_back(f); // left corner
//			if (i==n-3) data_.push_back(f); // right corner
		}

		assert(data_.size()==n-1);
//		data_[0].adjustCorner(SymmetryFactorType::CORNER_LEFT);
//		data_[n-1].adjustCorner(SymmetryFactorType::CORNER_RIGHT);
	}

	const SymmetryFactorType& operator()(size_t site) const
	{
		assert(site<data_.size());
		return data_[site];
	}

private:

	ProgramGlobals::Vector<SymmetryFactorType>::Type data_;

}; // SymmetryLocal

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_LOCAL_H

