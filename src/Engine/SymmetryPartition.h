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

#ifndef SYMMETRY_PARTITION_H
#define SYMMETRY_PARTITION_H

#include "ProgramGlobals.h"
#include "IoSimple.h"

namespace Mpspp {

class SymmetryPartition {

	class SymmetryPartitionOpaque {

	public:

		SymmetryPartitionOpaque(size_t length,size_t offset)
			: length_(length),offset_(offset)
		{}

		size_t size() const { return length_; }

		size_t offset() const { return offset_; }

	private:

		size_t length_;
		size_t offset_;
	};

public:

	typedef PsimagLite::IoSimple::In IoInputType;

	SymmetryPartition(const ProgramGlobals::Vector<size_t>::Type& data)
	: data_(data)
	{}

	size_t size() const
	{
		return data_.size();
	}

//	size_t offset() const
//	{
//		return offset_;
//	}

	SymmetryPartitionOpaque operator()(size_t i) const
	{
		return SymmetryPartitionOpaque(data_[i+1]-data_[i],data_[i]);
	}

private:

	const ProgramGlobals::Vector<size_t>::Type& data_;
//	size_t offset_;

}; // SymmetryPartition

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_PARTITION_H

