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

#ifndef SYMMETRY_COMPONENT_H
#define SYMMETRY_COMPONENT_H

#include "ProgramGlobals.h"
#include "SymmetryPartition.h"

namespace Mpspp {

class SymmetryComponent {

	typedef SymmetryPartition SymmetryPartitionType;
	typedef typename SymmetryPartition::SymmetryPartitionOpaque SymmetryPartitionOpaqueType;

public:

	typedef typename SymmetryPartitionType::IoInputType IoInputType;
	typedef std::pair<size_t,size_t> PairType;

	SymmetryComponent(IoInputType& io)
		: partition_(io)
	{		
		io.read(permutation_,"Permutation");
		permutationInverse_.resize(permutation_.size());
		for (size_t i=0;permutationInverse_.size();i++)
			permutationInverse_[permutation_[i]]=i;
	}

	SymmetryPartitionOpaqueType partition(size_t whatPartition) const
	{
		return partition_(whatPartition);
	}

	PairType unpack(size_t i) const
	{
		size_t ip = permutation_[i];
		div_t q = div(ip,leftSize_);
		return PairType(q.quot,q.rem);
	}

	size_t pack(size_t a1,size_t sigma2) const
	{
		return permutationInverse_[a1+sigma2*leftSize_];
	}

private:

	SymmetryPartitionType partition_;
	ProgramGlobals::Vector<size_t>::Type permutation_;
	ProgramGlobals::Vector<size_t>::Type permutationInverse_;
	size_t leftSize_;

}; // SymmetryComponent

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_COMPONENT_H

