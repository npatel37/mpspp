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
#include "IoSimple.h"

namespace Mpspp {

class SymmetryComponent {


public:

	typedef PsimagLite::IoSimple::In IoInputType;
	typedef std::pair<size_t,size_t> PairType;

	enum {CORNER_LEFT,CORNER_RIGHT};

	SymmetryComponent(IoInputType& io,size_t leftSize)
	: leftSize_(leftSize)
	{		
		loadInternal(io);
	}

	void adjustCorner(size_t corner)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to write adjustCorner. I cannot go further until this is implemented\n";
//		throw std::runtime_error(str.c_str());
		std::cerr<<str;
	}

	size_t partitionSize(size_t i) const
	{
		assert(i+1<partition_.size());
		return partition_[i+1]-partition_[i];
	}

	size_t partitionOffset(size_t i) const
	{
		assert(i<partition_.size());
		return partition_[i];
	}

	PairType unpack(size_t i) const
	{
		size_t ip = permutation_[i];
		div_t q = div(ip,leftSize_);
		return PairType(q.rem,q.quot);
	}

	size_t pack(size_t a1,size_t sigma2) const
	{
		return permutationInverse_[a1+sigma2*leftSize_];
	}

	size_t size() const
	{
		return quantumNumbers_.size();
	}

	const std::vector<size_t>& block() const { return block_; }

private:

	// Match with SaveInternal in DMRG++'s Basis.h
	template<typename IoInputter>
	void loadInternal(IoInputter& io)
	{
		//int x=0;
		//useSu2Symmetry_=false;
		//io.readline(x,"#useSu2Symmetry=");
		//if (x>0) useSu2Symmetry_=true;
		io.read(block_,"#BLOCK");
		io.read(quantumNumbers_,"#QN");
		//io.read(electrons_,"#ELECTRONS");
		//io.read(electronsOld_,"#0OLDELECTRONS");
		io.read(partition_,"#PARTITION");
		io.read(permutationInverse_,"#PERMUTATIONINVERSE");
		permutation_.resize(permutationInverse_.size());
		for (size_t i=0;i<permutation_.size();i++) permutation_[permutationInverse_[i]]=i;
		/*dmrgTransformed_=false;
		if (useSu2Symmetry_)
			symmSu2_.load(io);
		else
			symmLocal_.load(io);*/
	}

	std::vector<size_t> block_;
	std::vector<size_t> quantumNumbers_;
	std::vector<size_t> partition_;
	ProgramGlobals::Vector<size_t>::Type permutation_;
	ProgramGlobals::Vector<size_t>::Type permutationInverse_;
	size_t leftSize_;

}; // SymmetryComponent

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_COMPONENT_H

